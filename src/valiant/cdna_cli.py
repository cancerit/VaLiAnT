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
import os
import sys
from typing import Callable, Union, Dict, Tuple, Optional
import click
import numpy as np
import pandas as pd
from .constants import CDS_ONLY_MUTATORS
from .cli_utils import load_codon_table, validate_adaptor, set_logger
from .common_cli import common_params, existing_file
from .enums import TargetonMutator
from .models.base import StrandedPositionRange, PositionRange
from .models.cdna import CDNA, AnnotatedCDNA
from .models.cdna_seq_repository import CDNASequenceRepository
from .models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from .models.codon_table import CodonTable
from .models.mutated_sequences import MutationCollection
from .models.sequences import Sequence
from .models.snv_table import AuxiliaryTables
from .models.targeton import Targeton, CDSTargeton
from .utils import get_constant_category
from .writers import write_oligo_metadata


def get_cdna(
    seq_repo: CDNASequenceRepository,
    codon_table: CodonTable,
    targeton_cfg: CDNATargetonConfig
) -> CDNA:
    cdna_id: str = targeton_cfg.seq_id
    cdna = seq_repo.get(cdna_id)

    if not cdna:
        raise KeyError(f"Sequence identifier '{cdna_id}' not found!")

    if isinstance(cdna, AnnotatedCDNA):
        cdna.validate(codon_table)
    elif CDS_ONLY_MUTATORS & targeton_cfg.mutators:
        raise ValueError("Invalid mutator for non-annotated cDNA targeton!")

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
) -> Tuple[str, Dict[TargetonMutator, MutationCollection]]:
    sequence, cds_ext_5p, cds_ext_3p = cdna.get_extended_subsequence(
        targeton_cfg.r2_range)
    seq = sequence.sequence
    return seq, CDSTargeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range), cds_ext_5p, cds_ext_3p).compute_mutations(
                targeton_cfg.mutators, aux)


def get_cdna_mutations(
    aux: AuxiliaryTables,
    targeton_cfg: CDNATargetonConfig,
    cdna: AnnotatedCDNA
) -> Tuple[str, Dict[TargetonMutator, MutationCollection]]:
    seq = cdna.get_subsequence_string(targeton_cfg.r2_range)
    return seq, Targeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range)).compute_mutations(
                targeton_cfg.mutators)


def process_targeton(
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
    get_mutations = (
        get_annotated_cdna_mutations if isinstance(cdna, AnnotatedCDNA) else
        get_cdna_mutations
    )
    r2, mut_collections = get_mutations(aux, targeton_cfg, cdna)  # type: ignore

    # Get constant sequences (if any)
    c1 = cdna.get_subsequence_string(
        PositionRange(t_start, r2_start - 1)) if t_start != r2_start else ''
    c2 = cdna.get_subsequence_string(
        PositionRange(r2_end + 1, t_end)) if t_end != r2_end else ''

    # Merge mutation collections
    # TODO: add mutator
    df = pd.concat([mc.df for mc in mut_collections.values()])

    # Offset relative mutation position
    df.mut_position += 1

    df['ref_seq'] = get_constant_category(c1 + r2 + c2, df.shape[0])
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

    # TODO: apply maximum length mask
    df['oligo_length'] = df.mseq.str.len().astype(np.int32)

    return df


@click.command()
@common_params
@click.option('--annot', type=existing_file)
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
    targetons = CDNATargetonConfigCollection.load(oligo_info)

    # Load cDNA sequences
    seq_repo = CDNASequenceRepository.load(
        targetons.sequence_ids, ref_fasta, annot_fp=annot)

    # Get auxiliary tables
    aux: AuxiliaryTables = get_auxiliary_tables(
        targetons, load_codon_table(codon_table))

    get_cdna_f = partial(get_cdna, seq_repo, aux.codon_table)
    for targeton_cfg in targetons.cts:

        # Generate oligonucleotides
        meta = process_targeton(
            get_cdna_f, aux, adaptor_5, adaptor_3, targeton_cfg)

        # Set metadata
        meta['species'] = get_constant_category(species, meta.shape[0])
        meta['assembly'] = get_constant_category(assembly, meta.shape[0])

        # TODO: define file name format
        fn = f"{targeton_cfg.seq_id}_{hex(hash(targeton_cfg))}.csv"
        write_oligo_metadata(meta, os.path.join(output, fn))
