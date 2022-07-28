########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#############################

from dataclasses import dataclass
from functools import partial
import logging
import sys
from typing import Callable, Dict, Optional, NoReturn
import click
import numpy as np
import pandas as pd

from .cli_utils import load_codon_table, log_excluded_oligo_counts
from .common_cli import common_params, existing_file
from .constants import CDS_ONLY_MUTATORS, META_MAVE_NT, META_MSEQ, META_MSEQ_NO_ADAPT, META_OLIGO_NAME, MUTATOR_CATEGORIES
from .enums import TargetonMutator
from .errors import SequenceNotFound, InvalidMutatorForTarget
from .metadata_utils import get_mave_nt_from_row
from .models.base import StrandedPositionRange, PositionRange
from .models.cdna import CDNA, AnnotatedCDNA
from .models.cdna_seq_repository import CDNASequenceRepository
from .models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from .models.codon_table import CodonTable
from .models.cdna_config import CDNAConfig
from .models.metadata_table import MetadataTable
from .models.mutated_sequences import MutationCollection
from .models.oligo_generation_info import OligoGenerationInfo
from .models.oligo_renderer import get_oligo_names_from_dataframe
from .models.oligo_template import MUTATION_TYPE_CATEGORIES_T, decode_mut_types_cat
from .models.options import Options
from .models.snv_table import AuxiliaryTables
from .models.targeton import Targeton, CDSTargeton
from .utils import get_constant_category, get_empty_category_column, get_source_type_column, repr_enum_list


def exit_on_critical_exception(ex: Exception, msg: str) -> NoReturn:
    logging.critical(ex.args[0])
    logging.critical(msg)
    sys.exit(1)


def get_cdna(
    seq_repo: CDNASequenceRepository,
    targeton_cfg: CDNATargetonConfig
) -> CDNA:
    cdna_id: str = targeton_cfg.seq_id
    cdna = seq_repo.get(cdna_id)

    if not cdna:
        raise KeyError(f"Sequence identifier '{cdna_id}' not found!")

    cds_only_mutators = targeton_cfg.mutators & CDS_ONLY_MUTATORS
    if cds_only_mutators:
        if isinstance(cdna, AnnotatedCDNA):

            # CDS-only mutators in non-CDS target
            if targeton_cfg.r2_range not in cdna.cds_range:
                raise InvalidMutatorForTarget(
                    "Invalid mutators for non-CDS cDNA targeton: "
                    f"{repr_enum_list(cds_only_mutators)}!")

        else:

            # CDS-only mutators when the annotation is not available
            raise InvalidMutatorForTarget(
                "Invalid mutators for non-annotated cDNA targeton: "
                f"{repr_enum_list(cds_only_mutators)}!")

    return cdna


def get_auxiliary_tables(
    targetons: CDNATargetonConfigCollection,
    codon_table: CodonTable
) -> AuxiliaryTables:
    all_mutators: frozenset = targetons.mutators
    return AuxiliaryTables(
        codon_table,
        frozenset(['+']),
        TargetonMutator.SNV in all_mutators,
        TargetonMutator.SNV_RE in all_mutators,
        TargetonMutator.AA in all_mutators)


def get_annotated_cdna_mutations(
    aux: AuxiliaryTables,
    targeton_cfg: CDNATargetonConfig,
    cdna: AnnotatedCDNA
) -> Dict[TargetonMutator, MutationCollection]:
    sequence, cds_ext_5p, cds_ext_3p = cdna.get_extended_subsequence(
        targeton_cfg.r2_range)
    seq = sequence.sequence
    spr = StrandedPositionRange.to_plus_strand(targeton_cfg.r2_range)
    targeton = CDSTargeton.build_without_variants(spr, seq, cds_ext_5p, cds_ext_3p)
    return targeton.compute_mutations(targeton_cfg.mutators, aux)


def get_cdna_mutations(
    targeton_cfg: CDNATargetonConfig,
    cdna: CDNA
) -> Dict[TargetonMutator, MutationCollection]:
    seq = cdna.get_subsequence_string(targeton_cfg.r2_range)
    spr = StrandedPositionRange.to_plus_strand(targeton_cfg.r2_range)
    targeton = Targeton.build_without_variants(spr, seq)
    return targeton.compute_mutations(targeton_cfg.mutators)


def mut_coll_to_df(
    get_empty_aa_column: Callable[[int], pd.Categorical],
    r2_start: int,
    mutator: TargetonMutator,
    mc: MutationCollection
) -> pd.DataFrame:

    # Beware: the dataframe will be modified in place
    df = mc.df
    if df is None:
        raise ValueError("Missing dataframe!")

    rown = df.shape[0]

    # Offset relative mutation position
    df.mut_position += r2_start

    # Populate mutator field
    df['mutator'] = get_constant_category(
        mutator.value, rown, categories=MUTATOR_CATEGORIES)

    # Avoid categorical to object conversion on concatenation
    for aa_field in ['ref_aa', 'alt_aa']:
        if aa_field not in df.columns:
            df[aa_field] = get_empty_aa_column(rown)

    # Decode mutation type (if set)
    df['mut_type'] = (
        decode_mut_types_cat(df.mut_type) if 'mut_type' in df.columns else
        get_empty_category_column(MUTATION_TYPE_CATEGORIES_T, rown)
    )

    return df


def get_targeton_metadata_table(
    get_cdna_f: Callable,
    aux: AuxiliaryTables,
    adaptor_5: Optional[str],
    adaptor_3: Optional[str],
    targeton_cfg: CDNATargetonConfig
) -> pd.DataFrame:
    cdna = get_cdna_f(targeton_cfg)
    t_start, t_end = targeton_cfg.targeton_range.to_tuple()
    r2_start, r2_end = targeton_cfg.r2_range.to_tuple()

    # Generate mutated sequences
    mut_collections = (
        get_annotated_cdna_mutations(aux, targeton_cfg, cdna) if isinstance(cdna, AnnotatedCDNA) else
        get_cdna_mutations(targeton_cfg, cdna)
    )

    # Get constant sequences (if any)
    c1 = cdna.get_subsequence_string(
        PositionRange(t_start, r2_start - 1)) if t_start != r2_start else ''
    c2 = cdna.get_subsequence_string(
        PositionRange(r2_end + 1, t_end)) if t_end != r2_end else ''

    ref_seq = cdna.get_subsequence_string(targeton_cfg.targeton_range)

    get_empty_aa_column = partial(
        get_empty_category_column, tuple(aux.codon_table.amino_acid_symbols))

    # Merge mutation collections
    # Beware: the mutation collection dataframes will be modified
    df = pd.concat((
        mut_coll_to_df(get_empty_aa_column, r2_start, mutator, mc)
        for mutator, mc in mut_collections.items()
    ), ignore_index=True)
    del mut_collections

    rown: int = df.shape[0]

    # Compress reference and alternative sequences
    df.ref = df.ref.astype('category')
    df.new = df.new.astype('category')

    # Populate targeton attributes
    df['ref_seq'] = get_constant_category(ref_seq, df.shape[0])
    df['ref_start'] = np.int32(t_start)
    df['ref_end'] = np.int32(t_end)
    df['revc'] = np.int8(0)

    # Add constant sequences (if any)
    if c1:
        df[META_MSEQ] = c1 + df[META_MSEQ]
    if c2:
        df[META_MSEQ] = df[META_MSEQ] + c2

    df[META_MSEQ_NO_ADAPT] = df[META_MSEQ]

    # Add adaptors (if any)
    if adaptor_5:
        df[META_MSEQ] = adaptor_5 + df[META_MSEQ]

    if adaptor_3:
        df[META_MSEQ] = df[META_MSEQ] + adaptor_3

    # Add oligonucleotide lengths
    df['oligo_length'] = df.mseq.str.len().astype(np.int32)

    # Add sequence info (if available)
    if cdna.seq_info.gene_id:
        df['gene_id'] = get_constant_category(cdna.seq_info.gene_id, rown)
    if cdna.seq_info.transcript_id:
        df['transcript_id'] = get_constant_category(cdna.seq_info.transcript_id, rown)

    gene_label = cdna.seq_info.gene_id or 'NO_GENE'
    transcript_label = cdna.seq_info.transcript_id or 'NO_TRANSCRIPT'
    oligo_name_prefix = f"{targeton_cfg.seq_id}_{transcript_label}.{gene_label}_"

    # Generate oligonucleotide names
    df[META_OLIGO_NAME] = get_oligo_names_from_dataframe(oligo_name_prefix, df)

    # Add mutation MAVE-HGVS code
    df[META_MAVE_NT] = pd.Series(df.apply(get_mave_nt_from_row, axis=1), dtype='string')

    # Drop field that would be discarded downstream
    df = df.drop('var_type', axis=1)

    # Set sequence source type
    df['src_type'] = get_source_type_column('cdna', rown)

    return df


def process_targeton(
    species: str,
    assembly: str,
    aux: AuxiliaryTables,
    adaptor_5: Optional[str],
    adaptor_3: Optional[str],
    get_cdna_f: Callable,
    options: Options,
    output: str,
    targeton_cfg: CDNATargetonConfig
) -> OligoGenerationInfo:

    # Generate oligonucleotides
    metadata = MetadataTable.from_partial(
        species,
        assembly,
        get_targeton_metadata_table(
            get_cdna_f, aux, adaptor_5, adaptor_3, targeton_cfg),
        options.oligo_min_length,
        options.oligo_max_length)

    base_fn = f"{targeton_cfg.seq_id}_{targeton_cfg.get_hash()}"
    metadata.write_common_files(output, base_fn)

    # Log
    if metadata.short_oligo_n == 0:
        logging.warning(
            "Empty metadata table for cDNA targeton %s (%d-%d, %d-%d): no file generated!" % (
                targeton_cfg.seq_id,
                *targeton_cfg.targeton_range.to_tuple(),
                *targeton_cfg.r2_range.to_tuple()
            ))
        logging.warning(
            "Empty unique oligonucleotides table for cDNA targeton %s (%d-%d, %d-%d): no file generated!" % (
                targeton_cfg.seq_id,
                *targeton_cfg.targeton_range.to_tuple(),
                *targeton_cfg.r2_range.to_tuple()
            ))

    return metadata.get_info()


def run_cdna(config: CDNAConfig) -> None:
    options = config.get_options()

    # Load targeton configurations
    try:
        targetons = CDNATargetonConfigCollection.load(config.oligo_info_fp)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA targeton file!")

    # Load cDNA sequences
    try:
        seq_repo = CDNASequenceRepository.load(
            targetons.sequence_ids, config.ref_fasta_fp, annot_fp=config.annot_fp)
    except SequenceNotFound:
        sys.exit(1)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA annotation file!")

    # Get auxiliary tables
    aux: AuxiliaryTables = get_auxiliary_tables(
        targetons, load_codon_table(config.codon_table_fp))

    # Long oligonucleotides counter
    short_oligo_n: int = 0
    long_oligo_n: int = 0

    get_cdna_f = partial(get_cdna, seq_repo)
    process_targeton_f = partial(
        process_targeton,
        config.species,
        config.assembly,
        aux,
        config.adaptor_5,
        config.adaptor_3,
        get_cdna_f,
        options,
        config.output_dir)

    info: OligoGenerationInfo
    for targeton_cfg in targetons.cts:
        try:
            info = process_targeton_f(targeton_cfg)
            short_oligo_n += info.short_oligo_n
            long_oligo_n += info.long_oligo_n

        except ValueError as ex:
            exit_on_critical_exception(ex, "Failed to generate oligonucleotides!")

    log_excluded_oligo_counts(options, short_oligo_n, long_oligo_n)


@click.command()
@common_params
@click.option('--annot', type=existing_file, help="cDNA annotation file path")
def cdna(

    # Input files
    oligo_info: str,
    ref_fasta: str,
    codon_table: Optional[str],
    annot: Optional[str],

    # Output directory
    output: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: Optional[str],
    adaptor_3: Optional[str],

    # Actions
    max_length: int,
    min_length: int,

):
    """
    cDNA DMS oligonucleotide generation tool

    \b
    OLIGO_INFO is the BED-like input file path
    REF_FASTA is the cDNA multi-FASTA file path
    OUTPUT is the output directory path
    SPECIES will be included in the metadata
    ASSEMBLY will be included in the metadata
    """

    run_cdna(CDNAConfig(
        species=species,
        assembly=assembly,
        adaptor_5=adaptor_5,
        adaptor_3=adaptor_3,
        min_length=min_length,
        max_length=max_length,
        codon_table_fp=codon_table,
        oligo_info_fp=oligo_info,
        ref_fasta_fp=ref_fasta,
        output_dir=output,
        annot_fp=annot))
