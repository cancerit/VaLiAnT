########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from .variant import BaseVariantT, sort_variants


def compute_genomic_offset(variants: List[BaseVariantT]) -> int:
    """Given a list of variants, get the offset they would introduce if applied"""

    # TODO: check for overlapping variants
    return sum(
        variant.alt_ref_delta
        for variant in variants
    )


def _filter_variants_by_range(
    ref_start: int,
    ref_length: int,
    variants: List[BaseVariantT],
    sort: bool = False
) -> List[BaseVariantT]:
    a: List[BaseVariantT] = [
        variant
        for variant in variants
        if 0 <= (variant.genomic_position.position - ref_start) < ref_length
    ]
    return sort_variants(a) if sort else a


def _compute_alt_offsets(ref_start: int, ref_length: int, variants_in_range: List[BaseVariantT]) -> np.ndarray:
    """
    Generate an array of position offsets to convert relative positions in ALT to
    reference positions

    The variants are assumed to be in range and sorted by position.
    """

    alt_length: int = ref_length + compute_genomic_offset(variants_in_range)

    alt_offsets: np.ndarray = np.zeros(alt_length, dtype=np.int32)
    alt_offset: int = 0
    for variant in variants_in_range:

        # Filter out variants that do not add to ALT
        if variant.alt_ref_delta > 0:
            ref_offset: int = variant.genomic_position.position - ref_start

            # Increment offsets spanned by ALT
            for i in range(variant.alt_length):
                alt_offset += 1
                alt_offsets[ref_offset + i] += i

            # Set all downstream offsets (assumes the variants to be sorted by position)
            alt_offsets[ref_offset + variant.alt_length:] = alt_offset

    return alt_offsets


def _compute_ref_offsets(variants_in_range: List[BaseVariantT]) -> List[Tuple[int, int]]:
    """
    Map sorted positions to the cumulative offsets the corresponding variants introduce

    The variants are assumed to be in range and sorted by position.
    """

    offset: int = 0
    pos_offset = []
    for variant in variants_in_range:
        offset += variant.alt_ref_delta
        pos_offset.append((variant.genomic_position.position, offset))
    return pos_offset


@dataclass
class GenomicPositionOffsets:
    ref_start: int
    ref_length: int
    variants_in_range: List[BaseVariantT]
    _pos_offsets: List[Tuple[int, int]] = field(init=False)
    _ins_offsets: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self._pos_offsets = self._compute_ref_offsets()
        self._ins_offsets = self._compute_ins_offsets()

    def _compute_ins_offsets(self) -> List[Tuple[int, int]]:
        return _compute_alt_offsets(
            self.ref_start, self.ref_length, self.variants_in_range)

    def _compute_ref_offsets(self) -> List[Tuple[int, int]]:
        return _compute_ref_offsets(self.variants_in_range)

    @classmethod
    def from_variants(cls, ref_start: int, ref_length: int, variants: List[BaseVariantT]):
        variants_in_range = _filter_variants_by_range(
            ref_start, ref_length, variants, sort=True)
        return cls(ref_start, ref_length, variants_in_range)

    def get_offset(self, pos: int) -> int:
        if not self._pos_offsets or pos < self._pos_offsets[0]:
            return 0

        for p, offset in self._pos_offsets:
            if pos >= p:
                return offset
