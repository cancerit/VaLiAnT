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
from contextlib import nullcontext
import numpy as np
import pandas as pd
import pytest
from valiant.cdna_cli import get_cdna, get_auxiliary_tables, get_annotated_cdna_mutations, get_cdna_mutations, mut_coll_to_df, get_targeton_metadata_table
from valiant.enums import TargetonMutator
from valiant.models.base import PositionRange, StrandedPositionRange
from valiant.models.cdna import CDNA, AnnotatedCDNA
from valiant.models.cdna_seq_repository import CDNASequenceRepository
from valiant.models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from valiant.models.mutated_sequences import MutationCollection
from valiant.models.sequence_info import SequenceInfo
from valiant.models.sequences import Sequence
from valiant.models.snv_table import AuxiliaryTables
from valiant.utils import get_empty_category_column
from .utils import load_codon_table


POS_STRAND_SET = frozenset(['+'])

cdna_1_id = 'cdna_1'
cdna_1 = CDNA(Sequence("ATGACCGGTTTATACTGAGGTGTGTGTA"), SequenceInfo.empty())
cdna_1_annot = AnnotatedCDNA(cdna_1.seq, cdna_1.seq_info, StrandedPositionRange(1, 18, '+'))
cdna_id_map = {
    cdna_1_id: cdna_1
}
cdna_annot_id_map = {
    cdna_1_id: cdna_1_annot
}

targeton_metadata_table_fields = set([
    'mut_position',
    'ref',
    'new',
    'mseq',
    'mutator',
    'ref_aa',
    'alt_aa',
    'mut_type',
    'ref_seq',
    'ref_start',
    'ref_end',
    'revc',
    'oligo_length',
    'oligo_name',
    'src_type'
])

codon_table = load_codon_table()
aux = AuxiliaryTables(codon_table, POS_STRAND_SET, True, True, True)
get_empty_aa_column = partial(
    get_empty_category_column, tuple(aux.codon_table.amino_acid_symbols))


@pytest.mark.parametrize('is_annotated,pr,mutators,is_valid', [
    (False, cdna_1.range, frozenset([TargetonMutator.SNV]), True),
    (False, cdna_1.range, frozenset([TargetonMutator.SNV]), True),
    (False, StrandedPositionRange(19, 21, '+'), frozenset([TargetonMutator.SNV_RE]), False),  # Invalid mutator
    (True, StrandedPositionRange(19, 21, '+'), frozenset([TargetonMutator.SNV_RE]), False)  # Invalid mutator
])
def test_get_cdna(is_annotated, pr, mutators, is_valid):
    targeton_cfg = CDNATargetonConfig(cdna_1_id, cdna_1.range, pr, mutators)
    cdna_seq_repo = CDNASequenceRepository(
        cdna_annot_id_map if is_annotated else cdna_id_map)

    with pytest.raises(ValueError) if not is_valid else nullcontext():
        cdna = get_cdna(cdna_seq_repo, targeton_cfg)
        assert cdna == (cdna_1 if not is_annotated else cdna_1_annot)


def test_get_auxiliary_tables():
    pr = PositionRange(1, 10)
    aux = get_auxiliary_tables(CDNATargetonConfigCollection([
        CDNATargetonConfig('id', pr, pr, frozenset([
            TargetonMutator.SNV,
            TargetonMutator.SNV_RE,
            TargetonMutator.AA
        ]))
    ]), codon_table)

    # Validate SNV table
    assert aux._snv_table is not None
    assert aux._snv_table.strands == POS_STRAND_SET

    # Validate SNVRE table
    assert aux._snvre_table is not None
    assert aux._snvre_table._snv_table.strands == POS_STRAND_SET
    assert frozenset(aux._snvre_table._strand_all_syn_table.keys()) == POS_STRAND_SET
    assert frozenset(aux._snvre_table._strand_top_syn_table.keys()) == POS_STRAND_SET

    # Validate all amino acid table
    assert aux._all_aa_table is not None
    assert frozenset(aux._all_aa_table._strand_table.keys()) == POS_STRAND_SET


def test_get_annotated_cdna_mutations():
    pr = PositionRange(1, 10)
    mutators = frozenset([TargetonMutator.DEL1, TargetonMutator.SNV])
    targeton_cfg = CDNATargetonConfig('id', pr, pr, mutators)
    d = get_annotated_cdna_mutations(aux, targeton_cfg, cdna_1_annot)
    assert frozenset(d.keys()) == mutators
    for mutator in mutators:
        assert d[mutator].df is not None


def test_get_cdna_mutations():
    pr = PositionRange(1, 10)
    mutators = frozenset([TargetonMutator.DEL1, TargetonMutator.SNV])
    targeton_cfg = CDNATargetonConfig('id', pr, pr, mutators)
    d = get_cdna_mutations(targeton_cfg, cdna_1)
    assert frozenset(d.keys()) == mutators
    for mutator in mutators:
        assert d[mutator].df is not None


def test_mut_coll_to_df():
    raw_df = pd.DataFrame.from_records([
        (2, 0, 'A', 'T', 'TTGACCGGTT'),
        (2, 0, 'A', 'G', 'GTGACCGGTT'),
        (2, 0, 'A', 'C', 'CTGACCGGTT')
    ], columns=[
        'var_type',
        'mut_position',
        'ref',
        'new',
        'mseq'
    ])
    start_pos = 10

    df = mut_coll_to_df(
        get_empty_aa_column,
        start_pos,
        TargetonMutator.SNV,
        MutationCollection(df=raw_df.copy(deep=True)))

    # Check mutation position field
    assert df.mut_position.min() == raw_df.mut_position.min() + start_pos

    # Check amino acid fields
    assert 'ref_aa' in df
    assert 'alt_aa' in df

    # Check mutator field
    assert 'mutator' in df
    assert list(df.mutator.unique()) == [TargetonMutator.SNV.value]

    # Check mutation type field
    assert 'mut_type' in df


@pytest.mark.parametrize('annotated', [True, False])
def test_get_targeton_metadata_table(annotated):
    pr = PositionRange(1, 18)
    mutators = frozenset([TargetonMutator.SNV, TargetonMutator.DEL1])
    targeton_cfg = CDNATargetonConfig(cdna_1_id, pr, pr, mutators)

    adaptor_5 = 'TTTAAACCC'
    adaptor_3 = 'CCCAAATTT'
    repository = CDNASequenceRepository(
        cdna_annot_id_map if annotated else cdna_id_map)
    get_cdna_f = partial(get_cdna, repository)
    df = get_targeton_metadata_table(
        get_cdna_f, aux, adaptor_5, adaptor_3, targeton_cfg)

    assert set(df.columns) == targeton_metadata_table_fields

    ref_seqs = df.ref_seq.unique()
    assert len(ref_seqs) == 1
    assert ref_seqs[0] == repository.get_subsequence(cdna_1_id, pr).sequence

    # Check positions
    assert list(df.ref_start.unique()) == [pr.start]
    assert list(df.ref_end.unique()) == [pr.end]

    # Check oligonucleotide sequence
    assert df.mseq.str.startswith(adaptor_5).all()
    assert df.mseq.str.endswith(adaptor_3).all()

    # Check source type
    assert list(df.src_type.unique()) == ['cdna']

    # Check oligonucleotide length
    assert np.array_equal(df.oligo_length, df.mseq.str.len())
