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
    'src_type',
    'mave_nt'
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
        MutationCollection(raw_df.copy(deep=True)))

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
