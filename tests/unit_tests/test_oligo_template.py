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
from valiant.models.base import GenomicPosition, GenomicRange, TranscriptInfo
from valiant.models.new_pam import CdsPamBgAltSeqBuilder, PamBgAltSeqBuilder
from valiant.models.oligo_template import OligoTemplate, TargetonOligoSegment
from valiant.models.pam_protection import PamVariant
from valiant.models.refseq_ranges import ReferenceSequenceRanges
from valiant.models.sequences import ReferenceSequence
from valiant.models.targeton import PamProtCDSTargeton
from valiant.models.variant import CustomVariant, SubstitutionVariant
from .constants import DUMMY_PAM_PROTECTION_NT, DUMMY_PAM_SGRNA_ID
from .utils import get_aux_tables, get_pam_bg_alt_seq_builder, get_ref_seq

TRANSCRIPT_INFO = TranscriptInfo('GENE_ID', 'TRANSCRIPT_ID')


aux = get_aux_tables()


@pytest.mark.parametrize('targetons', [
    ['GGGACGTACGTTTTTACGTGGG', 'ATATATATATATATATATATATA', 'GGGGTATGCGGATTGACTGCTAGT']
])
@pytest.mark.parametrize('pam_protection', [True, False])
@pytest.mark.parametrize('mutator', [TargetonMutator.DEL1, TargetonMutator.SNV])
def test_oligo_compute_mutations(targetons, mutator, pam_protection):
    gr = GenomicRange('X', 1, sum(len(seq) for seq in targetons), '+')
    ref_seq = ReferenceSequence(''.join(targetons), gr)

    pam_ref_seq = get_pam_bg_alt_seq_builder(ref_seq, pam_protection)

    adaptor_5 = 'AAAAAA'
    adaptor_3 = 'AAAAAA'

    targeton_start = gr.start
    segments = []
    for targeton_seq in targetons:
        targeton_length = len(targeton_seq)
        segments.append(
            get_pam_bg_alt_seq_builder(
                get_ref_seq(
                    targeton_seq,
                    chromosome=gr.chromosome,
                    pos=targeton_start,
                    strand=gr.strand
                ), pam_protection))
        targeton_start += targeton_length

    empty = frozenset()
    rsr = ReferenceSequenceRanges(
        'X', '+', 1, 22, 1, 22, (0, 0), (empty, empty, empty), empty)

    ot = OligoTemplate(rsr, TRANSCRIPT_INFO, pam_ref_seq, frozenset({DUMMY_PAM_SGRNA_ID}), set(), adaptor_5, adaptor_3, segments)
    for _, target_segment in ot.target_segments:
        mutation_collections = target_segment.compute_mutations(aux)
        df = mutation_collections[mutator].df
        if df is None:
            raise ValueError("Null data frame!")

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
    pam_ref_seq = PamBgAltSeqBuilder.from_ref_seq(
        ReferenceSequence(ref_seq, gr), [], list(pam_variants))

    print(pam_ref_seq)

    targeton = PamProtCDSTargeton(CdsPamBgAltSeqBuilder.from_noncds(
        pam_ref_seq, '', ''))

    sgrna_ids = frozenset({x.sgrna_id for x in pam_variants})
    segment = TargetonOligoSegment(targeton, frozenset())
    ot = OligoTemplate(
        rsr, TRANSCRIPT_INFO, pam_ref_seq, sgrna_ids, {custom_variant}, None, None, [segment])

    cvm = ot._get_custom_variant_mutation(custom_variant)

    assert cvm.pam_ref == pam_ref
