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
from valiant.models.annotated_sequence import CDSAnnotatedSequencePair
from valiant.models.base import GenomicPosition, GenomicRange, TranscriptInfo
from valiant.models.codon_table import CodonTable
from valiant.models.oligo_template import OligoTemplate, TargetonOligoSegment
from valiant.models.pam_protection import PamProtectedReferenceSequence, PamVariant
from valiant.models.refseq_ranges import ReferenceSequenceRanges
from valiant.models.sequences import Sequence, ReferenceSequence
from valiant.models.targeton import PamProtCDSTargeton, Targeton
from valiant.models.variant import CustomVariant, SubstitutionVariant
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

    empty = frozenset()
    rsr = ReferenceSequenceRanges(
        'X', '+', 1, 22, 1, 22, (0, 0), (empty, empty, empty), empty)

    ot = OligoTemplate(rsr, TRANSCRIPT_INFO, pam_ref_seq, set(), set(), adaptor_5, adaptor_3, segments)
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


def test_oligo_template_get_custom_variant_mutation():
    empty = frozenset()
    rsr = ReferenceSequenceRanges(
        'X', '+', 1, 22, 1, 22, (0, 0), (empty, empty, empty), empty)

    ref = 'A'
    pam_ref = 'C'
    mut_alt = 'T'
    pos = GenomicPosition('X', 3)
    pam_variant = PamVariant(pos, ref, pam_ref, 'sgrna-a')
    variant = SubstitutionVariant(pos, ref, mut_alt)

    custom_variant = CustomVariant(variant, None, None)
    pam_variants = {pam_variant}

    ref_seq = ref * 6

    gr = GenomicRange('X', 1, len(ref_seq), '+')
    pam_ref_seq = PamProtectedReferenceSequence.from_reference_sequence(
        ReferenceSequence(ref_seq, gr), pam_variants)

    print(pam_ref_seq)

    targeton = PamProtCDSTargeton.from_pam_seq(
        pam_ref_seq, '', '')

    sgrna_ids = {x.sgrna_id for x in pam_variants}
    segment = TargetonOligoSegment(targeton, set())
    ot = OligoTemplate(
        rsr, TRANSCRIPT_INFO, pam_ref_seq, sgrna_ids, {custom_variant}, None, None, [segment])

    cvm = ot._get_custom_variant_mutation(custom_variant)

    assert cvm.pam_ref == pam_ref
