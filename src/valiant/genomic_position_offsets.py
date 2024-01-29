########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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
from typing import Generic

from .uint_range import UIntRange
from .variant import VariantT, sort_variants


@dataclass(slots=True)
class PosOffset:
    """Relative position and cumulative offset starting from it"""

    pos: int
    offset: int


def compute_genomic_offset(variants: list[VariantT]) -> int:
    """Given a list of variants, get the offset they would introduce if applied"""

    # TODO: check for overlapping variants (best done upstream, when filtering the full list via PyRanges)
    return sum(
        variant.alt_ref_delta
        for variant in variants
    )


def _filter_variants_by_range(
    ref_start: int,
    ref_length: int,
    variants: list[VariantT],
    sort: bool = False
) -> list[VariantT]:
    a: list[VariantT] = [
        variant
        for variant in variants
        if 0 <= (variant.pos - ref_start) < ref_length
    ]
    return sort_variants(a) if sort else a


def _get_alt_length(ref_length: int, variants_in_range: list[VariantT]) -> int:
    alt_length = ref_length + compute_genomic_offset(variants_in_range)
    if alt_length < 0:
        raise ValueError("Negative ALT length: verify the variants!")
    return alt_length


def _compute_alt_offsets(ref_start: int, alt_length: int, variants_in_range: list[VariantT]) -> list[int]:
    """
    Generate an array of position offsets to convert relative positions in ALT to
    reference positions

    The variants are assumed to be in range and sorted by position.
    """

    # TODO: use an array instead (for performance)?
    alt_offsets = [0] * alt_length
    alt_offset: int = 0
    for variant in variants_in_range:

        # Filter out variants that do not add to ALT
        if variant.alt_ref_delta > 0:
            ref_offset: int = variant.pos - ref_start

            # Increment offsets spanned by ALT
            for i in range(variant.alt_len):
                alt_offset += 1
                alt_offsets[ref_offset + i] += i

            # Set all downstream offsets (assumes the variants to be sorted by position)
            offset = ref_offset + variant.alt_len
            alt_offsets[offset:] = [alt_offset] * (len(alt_offsets) - offset)

    return alt_offsets


def _compute_ref_offsets(variants_in_range: list[VariantT]) -> list[PosOffset]:
    """
    Map sorted positions to the cumulative offsets the corresponding variants introduce

    The variants are assumed to be in range and sorted by position.
    """

    offset: int = 0
    pos_offset = []
    for variant in variants_in_range:
        offset += variant.alt_ref_delta
        pos_offset.append(PosOffset(variant.pos, offset))
    return pos_offset


@dataclass(slots=True)
class GenomicPositionOffsets(Generic[VariantT]):
    ref_start: int
    ref_length: int
    variants_in_range: list[VariantT]
    _alt_length: int = field(init=False)
    _pos_offsets: list[PosOffset] = field(init=False)
    _ins_offsets: list[int] = field(init=False)

    @property
    def alt_length(self) -> int:
        return self._alt_length

    def __post_init__(self) -> None:
        if self.ref_start < 0:
            raise ValueError(f"Invalid REF start {self.ref_start}!")
        if self.ref_length < 0:
            raise ValueError(f"Invalid REF length {self.ref_length}!")
        self._alt_length = _get_alt_length(self.ref_length, self.variants_in_range)
        self._pos_offsets = self._compute_ref_offsets()
        self._ins_offsets = self._compute_ins_offsets()

    def _compute_ins_offsets(self) -> list[int]:
        return _compute_alt_offsets(
            self.ref_start, self._alt_length, self.variants_in_range)

    def _compute_ref_offsets(self) -> list[PosOffset]:
        return _compute_ref_offsets(self.variants_in_range)

    @classmethod
    def from_variants(cls, ref_start: int, ref_length: int, variants: list[VariantT]):
        variants_in_range = _filter_variants_by_range(
            ref_start, ref_length, variants, sort=True)
        return cls(ref_start, ref_length, variants_in_range)

    def get_offset(self, pos: int) -> int:
        if not self._pos_offsets or pos < self._pos_offsets[0].pos:
            return 0

        # TODO: optimise search in sorted list
        for p in self._pos_offsets:
            if pos >= p.pos:
                return p.offset

        return 0

    def alt_to_ref_position(self, alt_pos: int) -> int:
        if not (self.ref_start <= alt_pos < self.ref_start + self._alt_length):
            raise ValueError(f"Invalid ALT position {alt_pos}: out of bounds!")

        return alt_pos - self._ins_offsets[alt_pos - self.ref_start]

    def ref_to_alt_position(self, ref_pos: int) -> int:
        return ref_pos + self.get_offset(ref_pos)

    def ref_to_alt_range(self, r: UIntRange) -> UIntRange:
        return UIntRange(
            self.ref_to_alt_position(r.start),
            self.ref_to_alt_position(r.end))

    def ref_to_alt_variant(self, variant: VariantT) -> VariantT:
        return variant.clone(pos=self.ref_to_alt_position(variant.pos))
