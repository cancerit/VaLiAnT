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

from valiant.models.base import GenomicRange, GenomicPosition
from valiant.models.new_pam import CdsPamBgAltSeqBuilder, PamBgAltSeqBuilder
from valiant.models.refseq_ranges import ReferenceSequenceRanges, TargetonConfig
from valiant.models.sequences import ReferenceSequence
from valiant.models.oligo_segment import InvariantOligoSegment
from valiant.models.custom_variants import CustomVariant
from valiant.models.oligo_template import OligoTemplate
from valiant.models.targeton import PamProtTargeton, PamProtCDSTargeton
from valiant.models.pam_protection import PamVariant
from valiant.models.variant import DeletionVariant, SubstitutionVariant, InsertionVariant


CHR = 'X'
STRAND = '+'
SEQ_START = 500
EXON_REGION_INDEX = 1
EXON_PAM_OFFSET = 15


@pytest.mark.parametrize('variant', [
    DeletionVariant(GenomicPosition(CHR, SEQ_START + 1), 'A' * 25),
    SubstitutionVariant(GenomicPosition(CHR, SEQ_START + 1), 'A' * 25, 'C' * 20),
    InsertionVariant(GenomicPosition(CHR, SEQ_START + EXON_PAM_OFFSET), 'C' * 20)
])
@pytest.mark.parametrize('non_cds_region', [True, False])
def test_oligo_template_get_custom_variant_sgrna_ids(variant, non_cds_region):
    """
    Generate three regions, each with one PAM variant.
    Make all or only the second region exonic, depending on `non_cds_region`.
    """

    all_sgrna_ids = [f"sgrna-{i}" for i in range(3)]
    pam_variant_offsets = [5, EXON_PAM_OFFSET, 25]
    seq_length = 30
    seq_start = SEQ_START
    pam_variants = [
        PamVariant(GenomicPosition(CHR, seq_start + offset), 'A', 'T', sgrna_id)
        for sgrna_id, offset in zip(all_sgrna_ids, pam_variant_offsets)
    ]
    ref_seq = ReferenceSequence(
        'A' * seq_length, GenomicRange(CHR, seq_start, seq_start + seq_length - 1, STRAND))

    pam_seq = PamBgAltSeqBuilder.from_ref_seq(
        ref_seq, [], pam_variants)

    custom_variant = CustomVariant(variant, None, None)

    # Mark all or just the middle region as exonic
    def get_targeton(region_index: int, region: PamBgAltSeqBuilder):
        return (
            PamProtTargeton(region)
            if non_cds_region and region_index != EXON_REGION_INDEX else
            PamProtCDSTargeton(CdsPamBgAltSeqBuilder.from_noncds(region, 'AA', ''))
        )

    regions = [
        pam_seq.get_pam_sub(GenomicRange(CHR, start, end, STRAND))
        for start, end in [
            (seq_start + i * 10, seq_start + (i + 1) * 10 - 1)
            for i in range(3)
        ]
    ]

    segments = [
        InvariantOligoSegment(get_targeton(i, region))
        for i, region in enumerate(regions)
    ]

    empty = frozenset()
    rsr = ReferenceSequenceRanges.from_config(
        TargetonConfig('X', '+', 1, seq_length, 1, seq_length, (0, 0), (empty, empty, empty), empty))

    ot = OligoTemplate(
        rsr, None, None, pam_seq, frozenset(all_sgrna_ids), {custom_variant}, None, None, segments)

    # Get sgRNA of the PAM variants in codons that are included in the custom mutation
    sgrna_ids = ot._get_custom_variant_sgrna_ids(custom_variant)

    sgrna_ids_exp = frozenset(
        all_sgrna_ids if not isinstance(variant, InsertionVariant) else
        [all_sgrna_ids[EXON_REGION_INDEX]]
    )

    assert sgrna_ids == sgrna_ids_exp
