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
import pandas as pd
from typing import Callable, Union, Dict
from .constants import CDS_ONLY_MUTATORS
from .enums import TargetonMutator
from .models.base import StrandedPositionRange
from .models.cdna import CDNA, AnnotatedCDNA
from .models.cdna_seq_repository import CDNASequenceRepository
from .models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from .models.codon_table import CodonTable
from .models.mutated_sequences import MutationCollection
from .models.sequences import Sequence
from .models.snv_table import AuxiliaryTables
from .models.targeton import Targeton, CDSTargeton
from .cli import _load_codon_table


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
) -> Dict[TargetonMutator, MutationCollection]:
    sequence, cds_ext_5p, cds_ext_3p = cdna.get_extended_subsequence(
        targeton_cfg.r2_range)
    seq = sequence.sequence
    return CDSTargeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range), cds_ext_5p, cds_ext_3p).compute_mutations(
                targeton_cfg.mutators, aux)


def get_cdna_mutations(
    aux: AuxiliaryTables,
    targeton_cfg: CDNATargetonConfig,
    cdna: AnnotatedCDNA
) -> Dict[TargetonMutator, MutationCollection]:
    seq = cdna.get_subsequence_string(targeton_cfg.r2_range)
    return Targeton(
        seq, seq, StrandedPositionRange.to_plus_strand(
            targeton_cfg.r2_range)).compute_mutations(
                targeton_cfg.mutators)


def main(targeton_fp: str, fasta_fp: str, annot_fp: str) -> None:

   # Load targeton configurations
    targetons = CDNATargetonConfigCollection.load(targeton_fp)

    # Load cDNA sequences
    seq_repo = CDNASequenceRepository.load(
        targetons.sequence_ids, fasta_fp, annot_fp=annot_fp)

    # Get auxiliary tables
    codon_table = _load_codon_table(None)
    aux: AuxiliaryTables = get_auxiliary_tables(targetons, codon_table)

    f = partial(get_cdna, seq_repo, codon_table)
    for targeton_cfg in targetons.cts:
        cdna = f(targeton_cfg)

        # Get target sequence
        r2 = cdna.get_subsequence_string(targeton_cfg.r2_range)

        # Generate mutated sequences
        get_mutations = (
            get_annotated_cdna_mutations if isinstance(cdna, AnnotatedCDNA) else
            get_cdna_mutations
        )
        mut_collections = get_mutations(aux, targeton_cfg, cdna)  # type: ignore

        # Get constant sequences
        c1 = cdna.get_subsequence_string_before(targeton_cfg.r2_range.start)
        c2 = cdna.get_subsequence_string_after(targeton_cfg.r2_range.end)

        df = pd.concat([mc.df for mc in mut_collections.values()])
        df.mut_position += 1
        print(df)


if __name__ == '__main__':
    _, targeton_fp, fasta_fp, annot_fp = sys.argv
    main(targeton_fp, fasta_fp, annot_fp)
