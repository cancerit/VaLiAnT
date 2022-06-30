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

import pytest
from valiant.enums import TargetonMutator
from valiant.models.base import GenomicRange, TranscriptInfo
from valiant.models.codon_table import CodonTable
from valiant.models.oligo_template import OligoTemplate, TargetonOligoSegment
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.sequences import Sequence, ReferenceSequence
from valiant.models.targeton import Targeton
from .constants import CODON_TABLE_FP, DUMMY_PAM_PROTECTION_NT
from .utils import get_data_file_path, get_no_op_pam_protected_sequence, get_targeton

TRANSCRIPT_INFO = TranscriptInfo('GENE_ID', 'TRANSCRIPT_ID')


@pytest.mark.parametrize('targetons', [
    ['GGGACGTACGTTTTTACGTGGG', 'ATATATATATATATATATATATA', 'GGGGTATGCGGATTGACTGCTAGT']
])
@pytest.mark.parametrize('pam_protection', [True, False])
@pytest.mark.parametrize('mutator', [TargetonMutator.DEL1, TargetonMutator.SNV])
def test_oligo_compute_mutations(targetons, mutator, pam_protection):

    def get_segment(seq):
        return TargetonOligoSegment(
            get_targeton(seq, pam_protection), {mutator})

    ct = CodonTable.load(get_data_file_path(CODON_TABLE_FP))
    ref_seq = ''.join(targetons)

    gr = GenomicRange('X', 1, sum(len(seq) for seq in targetons), '+')
    pam_ref_seq = get_no_op_pam_protected_sequence(ref_seq, gr)

    adaptor_5 = 'AAAAAA'
    adaptor_3 = 'AAAAAA'
    segments = list(map(get_segment, targetons))

    ot = OligoTemplate(TRANSCRIPT_INFO, pam_ref_seq, set(), set(), adaptor_5, adaptor_3, segments)
    for _, target_segment in ot.target_segments:
        mutation_collections = target_segment.compute_mutations(ct)
        df = mutation_collections[mutator].df
        if pam_protection:
            if mutator == TargetonMutator.DEL1:
                assert all(
                    set(m.mseq) == {DUMMY_PAM_PROTECTION_NT}
                    for m in df.itertuples()
                )
            elif mutator == TargetonMutator.SNV:
                assert all(
                    set(m.mseq) - {m.new} == {DUMMY_PAM_PROTECTION_NT}
                    for m in df.itertuples()
                )
