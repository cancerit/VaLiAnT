########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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
from typing import Iterable, List, Optional
from valiant.models.background_variants import GenomicPositionOffsets
from valiant.models.base import GenomicRange
from valiant.models.new_pam import PamBgAltSeqBuilder
from valiant.models.pam_protection import PamVariant

from valiant.models.sequences import ReferenceSequence
from valiant.models.variant import BaseVariant


@dataclass(frozen=True)
class PamProtectedReferenceSequence(ReferenceSequence):
    __slots__ = [
        'sequence',
        'genomic_range',
        'alt_seq',
        'bg_offsets'
    ]

    alt_seq: PamBgAltSeqBuilder
    bg_offsets: GenomicPositionOffsets

    @property
    def pam_protected_sequence(self) -> str:
        # TODO: CDS extension?
        return self.alt_seq.get_pam_seq()

    @property
    def bg_sequence(self) -> str:
        return self.alt_seq.get_bg_seq()

    @property
    def pam_variants(self) -> List[PamVariant]:
        return self.alt_seq.pam_variants

    @property
    def bg_variants(self) -> List[BaseVariant]:
        return self.alt_seq.bg_variants

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.pam_protected_sequence) != len(self.genomic_range):
            raise ValueError("PAM protected sequence and genomic range have different lengths!")

    @classmethod
    def from_reference_sequence(
        cls,
        ref_seq: ReferenceSequence,
        pam_variants: Iterable[PamVariant],
        bg_variants: Iterable[BaseVariant] = list()
    ) -> 'PamProtectedReferenceSequence':
        alt_seq = PamBgAltSeqBuilder.from_ref_seq(ref_seq, bg_variants, pam_variants)

        # TODO: apply masking and link to region boundary protection
        gpo = GenomicPositionOffsets.from_variants(
            ref_seq.genomic_range.start, len(ref_seq), alt_seq.bg_variants)

        return cls(
            ref_seq.sequence,
            ref_seq.genomic_range,
            alt_seq,
            gpo)

    def get_subsequence(self, genomic_range: GenomicRange) -> 'PamProtectedReferenceSequence':
        ref_seq: ReferenceSequence = super().get_subsequence(genomic_range)
        pam_variants = [
            variant
            for variant in self.pam_variants
            if genomic_range.contains_position(variant.genomic_position)
        ]
        bg_variants = [
            variant
            for variant in self.bg_variants
            if genomic_range.contains_position(variant.genomic_position)
        ]
        return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_variants, bg_variants)

    def _validate_variant_position(self, variant: BaseVariant) -> None:
        if not self.genomic_range.contains_position(variant.genomic_position):
            raise ValueError("Variant not in genomic range!")

    def apply_variant(self, variant: BaseVariant, ref_check: bool = False) -> str:
        self._validate_variant_position(variant)
        return variant.mutate(self.pam_protected_sequence, self.genomic_range.start, ref_check=ref_check)

    def get_variant_corrected_ref(self, variant: BaseVariant) -> Optional[str]:
        self._validate_variant_position(variant)
        return variant.get_corrected_ref(self.pam_protected_sequence, self.genomic_range.start)
