########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
#
# The Software remains the property of Genome Research Ltd ("GRL").
#
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
#
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
#
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
#
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
#############################

from functools import partial
import logging
import sys
from typing import Callable, Dict, Optional, NoReturn
import click
import numpy as np
import pandas as pd
from .constants import CDS_ONLY_MUTATORS, MUTATOR_CATEGORIES
from .cli_utils import load_codon_table, validate_adaptor, set_logger
from .common_cli import common_params, existing_file
from .enums import TargetonMutator
from .errors import SequenceNotFound, InvalidMutatorForTarget
from .models.base import StrandedPositionRange, PositionRange
from .models.cdna import CDNA, AnnotatedCDNA
from .models.cdna_seq_repository import CDNASequenceRepository
from .models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from .models.codon_table import CodonTable
from .models.metadata_table import MetadataTable
from .models.mutated_sequences import MutationCollection
from .models.oligo_generation_info import OligoGenerationInfo
from .models.oligo_renderer import get_oligo_name
from .models.oligo_template import MUTATION_TYPE_CATEGORIES_T, decode_mut_types_cat
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
    return CDSTargeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range), cds_ext_5p, cds_ext_3p).compute_mutations(
                targeton_cfg.mutators, aux)


def get_cdna_mutations(
    targeton_cfg: CDNATargetonConfig,
    cdna: CDNA
) -> Dict[TargetonMutator, MutationCollection]:
    seq = cdna.get_subsequence_string(targeton_cfg.r2_range)
    return Targeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range)).compute_mutations(
                targeton_cfg.mutators)


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
    prefix = (adaptor_5 or '') + c1
    if prefix:
        df.mseq = prefix + df.mseq
    suffix = c2 + (adaptor_3 or '')
    if suffix:
        df.mseq = df.mseq + suffix

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
    df['oligo_name'] = pd.Series(
        df.apply(lambda r: get_oligo_name(
            oligo_name_prefix,
            r.var_type,
            r.mutator,
            r.mut_position,
            r.ref if not pd.isnull(r.ref) else None,
            r.new if not pd.isnull(r.new) else None), axis=1),
        dtype='string')

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
    max_length: int,
    output: str,
    targeton_cfg: CDNATargetonConfig
) -> OligoGenerationInfo:

    # Generate oligonucleotides
    metadata = MetadataTable.from_partial(
        species,
        assembly,
        get_targeton_metadata_table(
            get_cdna_f, aux, adaptor_5, adaptor_3, targeton_cfg),
        max_length)

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

    # Extra
    log: str

):
    # Set logging up
    set_logger(log)

    # Validate adaptor sequences
    validate_adaptor(adaptor_5)
    validate_adaptor(adaptor_3)

    # Load targeton configurations
    try:
        targetons = CDNATargetonConfigCollection.load(oligo_info)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA targeton file!")

    # Load cDNA sequences
    try:
        seq_repo = CDNASequenceRepository.load(
            targetons.sequence_ids, ref_fasta, annot_fp=annot)
    except SequenceNotFound:
        sys.exit(1)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA annotation file!")

    # Get auxiliary tables
    aux: AuxiliaryTables = get_auxiliary_tables(
        targetons, load_codon_table(codon_table))

    # Long oligonucleotides counter
    long_oligo_n: int = 0

    get_cdna_f = partial(get_cdna, seq_repo)
    process_targeton_f = partial(
        process_targeton,
        species,
        assembly,
        aux,
        adaptor_5,
        adaptor_3,
        get_cdna_f,
        max_length,
        output)

    for targeton_cfg in targetons.cts:
        try:
            long_oligo_n += process_targeton_f(targeton_cfg).long_oligo_n

        except ValueError as ex:
            exit_on_critical_exception(ex, "Failed to generate oligonucleotides!")

    # Log number of oligonucleotides discarded due to excessive length
    if long_oligo_n:
        logging.warning(
            "%d oligonucleotides longer than %d bases were discarded!" %
            (long_oligo_n, max_length))
