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

from array import array
from dataclasses import dataclass, field
from typing import Callable, Generic, Iterable

from .uint_range import UIntRange
from .utils import get_end
from .variant import VariantT, sort_variants
from .variant_group import VariantGroup


def _get_zeros(t: str, length: int) -> array:
    k = array(t).itemsize
    return array(t, bytes(k * length))


def get_u8_array(n: int) -> array:
    return _get_zeros('B', n)


def get_u32_array(n: int) -> array:
    return _get_zeros('I', n)


@dataclass(slots=True)
class PosOffset:
    """Relative position and cumulative offset starting from it"""

    pos: int
    offset: int


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


def _compute_alt_offsets(ref_start: int, alt_length: int, alt_variants: Iterable[VariantT]) -> tuple[array, array]:
    """
    Generate an array of position offsets to convert relative positions in ALT to
    reference positions

    The variants are assumed to be in range and sorted by position.
    """

    alt_offsets = get_u32_array(alt_length)
    ins_mask = get_u8_array(alt_length)
    alt_offset: int = 0
    for variant in alt_variants:

        # Filter out variants that do not add to ALT
        if variant.alt_ref_delta > 0:
            ref_offset: int = variant.pos - ref_start

            # Increment offsets spanned by ALT
            for i in range(variant.alt_len):
                offset = ref_offset + i
                alt_offset += 1
                alt_offsets[offset] += i
                ins_mask[offset] = 1

            # Set all downstream offsets (assumes the variants to be sorted by position)
            offset = ref_offset + variant.alt_len
            for i in range(offset, alt_length):
                alt_offsets[i] = alt_offset

    return alt_offsets, ins_mask


def _compute_ref_offsets(ref_variants: list[VariantT]) -> list[PosOffset]:
    """
    Map sorted positions to the cumulative offsets the corresponding variants introduce

    The variants are assumed to be in range and sorted by position.
    """

    offset: int = 0
    pos_offset = []
    for variant in ref_variants:
        offset += variant.alt_ref_delta
        pos_offset.append(PosOffset(variant.pos, offset))
    return pos_offset


def _build_variant_mask(origin: int, length: int, variants: Iterable[VariantT], cond: Callable[[VariantT], bool]) -> array:
    mask = get_u8_array(length)

    for variant in variants:

        # Filter out variants that do not add to ALT
        if cond(variant):
            ref_offset = variant.pos - origin

            # Increment offsets spanned by ALT
            for i in range(variant.alt_len):
                mask[ref_offset + i] = 1

    return mask


def _compute_ref_del_mask(ref_start: int, ref_length: int, ref_variants: list[VariantT]) -> array:
    return _build_variant_mask(
        ref_start, ref_length, ref_variants, lambda x: x.alt_ref_delta < 0)


@dataclass(slots=True)
class GenomicPositionOffsets(Generic[VariantT]):
    ref_start: int
    ref_length: int
    variants_in_range: VariantGroup[VariantT]

    # ALT
    _alt_length: int = field(init=False)

    # REF -> ALT
    _pos_offsets: list[PosOffset] = field(init=False)
    _ref_del_mask: array = field(init=False)

    # ALT -> REF
    _ins_offsets: array = field(init=False)
    _alt_ins_mask: array = field(init=False)

    @property
    def alt_length(self) -> int:
        return self._alt_length

    @property
    def alt_end(self) -> int | None:
        if self.alt_length == 0:
            return None

        return get_end(self.ref_start, self.alt_length)

    @property
    def alt_range(self) -> UIntRange | None:
        if self.alt_length == 0:
            return None
        return UIntRange.from_length(self.ref_start, self.alt_length)

    def __post_init__(self) -> None:
        if self.ref_start < 0:
            raise ValueError(f"Invalid REF start {self.ref_start}!")
        if self.ref_length < 0:
            raise ValueError(f"Invalid REF length {self.ref_length}!")

        # Compute the length of the ALT sequence
        self._alt_length = self.variants_in_range.get_alt_length(self.ref_length)
        if self._alt_length < 0:
            raise ValueError("Negative ALT length: verify the variants!")

        self._ref_del_mask = self._compute_del_mask()
        self._pos_offsets = self._compute_ref_offsets()
        self._ins_offsets, self._alt_ins_mask = self._compute_ins_offsets()

    def _compute_del_mask(self) -> array:
        return _compute_ref_del_mask(
            self.ref_start, self.ref_length, self.variants_in_range.variants)

    def _compute_ins_offsets(self) -> tuple[array, array]:
        return _compute_alt_offsets(
            self.ref_start, self._alt_length, self.variants_in_range.iter_alt_variants())

    def _compute_ref_offsets(self) -> list[PosOffset]:
        return _compute_ref_offsets(self.variants_in_range.variants)

    @classmethod
    def from_variants(cls, ref_start: int, ref_length: int, variants: list[VariantT]):
        variants_in_range = VariantGroup.from_variants(_filter_variants_by_range(
            ref_start, ref_length, variants, sort=True))
        return cls(ref_start, ref_length, variants_in_range)

    def get_offset(self, pos: int) -> int:
        if not self._pos_offsets or pos < self._pos_offsets[0].pos:
            return 0

        # TODO: optimise search in sorted list
        for p in self._pos_offsets:
            if pos >= p.pos:
                return p.offset

        return 0

    def validate_alt_position(self, alt_pos: int) -> None:
        if alt_pos < self.ref_start or self.alt_end is None or alt_pos > self.alt_end:
            raise ValueError(f"Invalid ALT position {alt_pos}: out of bounds!")

    def validate_alt_range(self, r: UIntRange) -> None:
        t = self.alt_range
        if not t or r not in t:
            raise ValueError(f"Invalid ALT range {r}: out of bounds!")

    def _pos_to_offset(self, pos: int) -> int:
        return pos - self.ref_start

    def alt_pos_exists_in_ref(self, alt_pos: int) -> bool:
        return self._alt_ins_mask[self._pos_to_offset(alt_pos)] == 0

    def ref_pos_exists_in_alt(self, ref_pos: int) -> bool:
        return self._ref_del_mask[self._pos_to_offset(ref_pos)] == 0

    def _alt_to_ref_position(self, alt_pos: int) -> int:
        """
        Lift a position from ALT to REF

        Unsafe, no boundary or insertion mask test.
        """

        return alt_pos - self._ins_offsets[self._pos_to_offset(alt_pos)]

    def alt_to_ref_position(self, alt_pos: int) -> int | None:
        """Lift a position from ALT to REF"""

        self.validate_alt_position(alt_pos)

        if not self.alt_pos_exists_in_ref(alt_pos):
            return None

        return self._alt_to_ref_position(alt_pos)

    def alt_range_to_ref_positions(self, r: UIntRange) -> list[int]:
        """Lift a range of positions from ALT to a sequence of positions in REF"""

        self.validate_alt_range(r)

        return [
            self._alt_to_ref_position(alt_pos)
            for alt_pos in range(r.start, r.start + len(r))
            if self.alt_pos_exists_in_ref(alt_pos)
        ]

    def ref_to_alt_position(self, ref_pos: int) -> int | None:
        if not self.ref_pos_exists_in_alt(ref_pos):
            return None
        return ref_pos + self.get_offset(ref_pos)

    def ref_to_alt_range(self, r: UIntRange) -> UIntRange:
        # TODO: consider whether there's a requirement for this
        #  to be able to return a range with gaps, considering
        #  that custom variants can't overlap background variants.
        start = self.ref_to_alt_position(r.start)
        end = self.ref_to_alt_position(r.end)
        assert start is not None
        assert end is not None
        return UIntRange(start, end)

    def ref_to_alt_variant(self, variant: VariantT) -> VariantT:
        return variant.clone(pos=self.ref_to_alt_position(variant.pos))
