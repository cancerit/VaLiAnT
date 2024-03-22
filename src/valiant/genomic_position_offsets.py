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
from dataclasses import dataclass
from enum import IntEnum
from typing import Callable, Iterable

from .array_utils import get_prev_index, get_next_index, get_u32_array, get_u8_array, get_i32_array
from .uint_range import UIntRange
from .utils import get_end
from .var_stats import VarStats, clamp_var_stats_collection, get_alt_ref_delta
from .variant import Variant, VariantT


class InvalidRefRange(ValueError):
    def __init__(self, r: UIntRange, *args: object) -> None:
        self.msg = f"Range {r} could not be lift over!"
        super().__init__(*args)


class SearchType(IntEnum):
    BEFORE = 0
    AFTER = 1


SEARCH_F: dict[SearchType, Callable[[array, int, int], int | None]] = {
    SearchType.BEFORE: get_prev_index,
    SearchType.AFTER: get_next_index
}


@dataclass(slots=True)
class PosOffset:
    """Relative position and cumulative offset starting from it"""

    pos: int
    offset: int


def get_pos_offset(pos_offsets: list[PosOffset], pos: int) -> int:
    if not pos_offsets:
        return 0

    first = pos_offsets[0]
    if pos < first.pos:
        return 0

    last = pos_offsets[-1]
    if pos >= last.pos:
        return last.offset

    for i in range(1, len(pos_offsets)):
        if pos < pos_offsets[i].pos:
            return pos_offsets[i - 1].offset

    return 0


def _compute_alt_ins_mask(ref_start: int, alt_length: int, alt_variants: Iterable[VarStats]) -> array:
    """
    Generate an array of position offsets to convert relative positions in ALT to
    reference positions

    The variants are assumed to be in range and sorted by position.
    """

    ins_mask = get_u8_array(alt_length)
    alt_offset: int = 0
    var_offset: int = 0
    for variant in alt_variants:

        # Filter out variants that do not add to ALT
        if variant.alt_ref_delta > 0:
            ref_offset: int = variant.pos + var_offset - ref_start

            alt_offset += var_offset

            # Increment offsets spanned by ALT
            for i in range(variant.alt_len):
                offset = ref_offset + i
                alt_offset += 1
                ins_mask[offset] = 1

            # Set all downstream offsets (assumes the variants to be sorted by position)
            offset = ref_offset + variant.alt_len

        var_offset += variant.alt_ref_delta

    return ins_mask


def _compute_ref_offsets(ref_variants: list[VarStats]) -> tuple[list[PosOffset], list[PosOffset]]:
    """
    Map sorted positions to the cumulative offsets the corresponding variants introduce

    The variants are assumed to be in range and sorted by position.
    """

    offset: int = 0
    pos_offset = []
    alt_offset = []
    alt_pos: int = 0
    for variant in ref_variants:
        if variant.alt_ref_delta != 0:
            alt_pos = variant.pos + offset
            offset += variant.alt_ref_delta
            pos_offset.append(PosOffset(variant.pos, offset))
            alt_offset.append(PosOffset(alt_pos, -offset))
    return pos_offset, alt_offset


def _compute_ref_del_mask(ref_start: int, ref_length: int, ref_variants: list[VarStats]) -> tuple[array, array]:
    del_mask = get_u8_array(ref_length)
    shift_mask = get_u8_array(ref_length)

    for variant in ref_variants:
        if variant.alt_ref_delta != 0:
            ref_offset = variant.pos - ref_start

            for i in range(max(1, variant.ref_len)):
                shift_mask[ref_offset + i] = 1

            # Filter out variants that do not add to ALT
            if variant.alt_ref_delta < 0:

                # Increment offsets spanned by ALT
                for i in range(variant.ref_len):
                    del_mask[ref_offset + i] = 1

    return del_mask, shift_mask


@dataclass(slots=True)
class GenomicPositionOffsets:
    ref_range: UIntRange

    # ALT
    alt_length: int

    # REF -> ALT
    _pos_offsets: list[PosOffset]
    _ref_del_mask: array
    _shift_mask: array

    # ALT -> REF
    _alt_offsets: list[PosOffset]
    _alt_ins_mask: array

    def _get_ref_pos_offset(self, ref_pos: int) -> int:
        return get_pos_offset(self._pos_offsets, ref_pos)

    def _get_alt_pos_offset(self, alt_pos: int) -> int:
        return get_pos_offset(self._alt_offsets, alt_pos)

    @classmethod
    def from_var_stats(cls, vs: Iterable[VarStats], r: UIntRange):
        cvs = clamp_var_stats_collection(vs, r)
        ref_start = r.start
        ref_length = len(r)

        alt_length = ref_length + get_alt_ref_delta(cvs)
        del_mask, shift_mask = _compute_ref_del_mask(ref_start, ref_length, cvs)
        pos_offsets, alt_offsets = _compute_ref_offsets(cvs)

        # TODO: build these based on the reference offsets
        ins_mask = _compute_alt_ins_mask(ref_start, alt_length, cvs)

        return cls(r, alt_length, pos_offsets, del_mask, shift_mask, alt_offsets, ins_mask)

    def __post_init__(self) -> None:
        if self.alt_length < 0:
            raise ValueError("Negative ALT length: verify the variants!")

    @property
    def ref_start(self) -> int:
        return self.ref_range.start

    @property
    def ref_length(self) -> int:
        return len(self.ref_range)

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

    def _ref_to_alt_offset(self, ref_pos: int) -> int:
        # TODO: optimise search in sorted list
        return self._get_ref_pos_offset(ref_pos)

    def get_offset(self, pos: int) -> int:
        if not self._pos_offsets or pos < self._pos_offsets[0].pos:
            return 0

        return self._ref_to_alt_offset(pos)

    def validate_ref_position(self, ref_pos: int) -> None:
        if ref_pos not in self.ref_range:
            raise ValueError(f"Invalid ALT position {ref_pos}: out of bounds!")

    def validate_alt_position(self, alt_pos: int) -> None:
        if alt_pos < self.ref_start or self.alt_end is None or alt_pos > self.alt_end:
            raise ValueError(f"Invalid ALT position {alt_pos}: out of bounds!")

    def validate_alt_range(self, r: UIntRange) -> None:
        t = self.alt_range
        if not t or r not in t:
            raise ValueError(f"Invalid ALT range {r}: out of bounds!")

    def _pos_to_offset(self, pos: int) -> int:
        return pos - self.ref_start

    def _offset_to_pos(self, offset: int) -> int:
        return self.ref_start + offset

    def alt_pos_exists_in_ref(self, alt_pos: int) -> bool:
        return self._alt_ins_mask[self._pos_to_offset(alt_pos)] == 0

    def ref_pos_exists_in_alt(self, ref_pos: int) -> bool:
        return self._ref_del_mask[self._pos_to_offset(ref_pos)] == 0

    def ref_pos_overlaps_var(self, ref_pos: int) -> bool:
        return self._shift_mask[self._pos_to_offset(ref_pos)] == 1

    def ref_var_overlaps_var(self, variant: Variant) -> bool:
        return variant.any_pos(self.ref_pos_overlaps_var)

    def alt_var_overlaps_var(self, variant: Variant) -> bool:
        # Lift start over
        ref_start = self.alt_to_ref_position(variant.pos)
        if ref_start is None:
            return True

        if variant.ref_len <= 1:
            return False

        # Lift end over
        ref_end = self.alt_to_ref_position(variant.ref_end)
        if ref_end is None:
            return True

        ref_range = UIntRange(ref_start, ref_end)

        # Test for change in reference length
        if len(ref_range) != variant.ref_len:
            return True

        # Test for lifted over positions overlapping indels
        for ref_pos in ref_range.positions:
            if self.ref_pos_overlaps_var(ref_pos):
                return True

        return False

    def _alt_to_ref_position(self, alt_pos: int) -> int:
        """
        Lift a position from ALT to REF

        Unsafe, no boundary or insertion mask test.
        """

        return alt_pos + self._get_alt_pos_offset(alt_pos)

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

    def _ref_offset_to_alt_pos(self, ref_pos: int) -> int:
        return ref_pos + self._ref_to_alt_offset(ref_pos)

    def ref_to_alt_position(self, ref_pos: int, nearest: SearchType | None = None) -> int | None:
        i = self._pos_to_offset(ref_pos)

        # Reference position before the context
        if i < 0:

            # This position is guaranteed to exist because alt positions are anchored to the left
            # Assumption: nonzero context length
            return ref_pos

        mask = self._ref_del_mask

        if i < self.ref_length:

            # Reference position within the context
            if mask[i] == 0:

                # A corresponding position exists
                return self._ref_offset_to_alt_pos(ref_pos)

            if nearest is None:
                return None

            # Search for an existing position before or after the query
            f = SEARCH_F[nearest]
            offset = f(mask, i, 0)
            if offset is None:
                return None

            # Return the nearest existing ALT position
            nearest_ref_pos = self._offset_to_pos(offset)
            return self._ref_offset_to_alt_pos(nearest_ref_pos)

        else:

            # Reference position after the context (guaranteed to exist)
            # Assumption: reference position is within the contig bounds (not verified anywhere)
            return self._ref_offset_to_alt_pos(self.ref_range.end)

    def ref_to_alt_range(self, r: UIntRange, shrink: bool = False) -> UIntRange | None:
        """
        Lift a range of positions from REF to ALT

        If either range boundary does not exist in ALT, and shrink is set,
        return the widest existing range within it.
        """

        # TODO: consider whether there's a requirement for this
        #  to be able to return a range with gaps, considering
        #  that custom variants can't overlap background variants.
        start = self.ref_to_alt_position(
            r.start, nearest=SearchType.AFTER if shrink else None)

        if start is None:
            return None

        end = self.ref_to_alt_position(
            r.end, nearest=SearchType.BEFORE if shrink else None)

        if end is None:
            raise RuntimeError("Unexpected null range end position!")

        if end < start:

            # No position in the range could be lift over
            return None

        return UIntRange(start, end)

    def ref_to_alt_variant(self, variant: VariantT) -> VariantT:
        return variant.clone(pos=self.ref_to_alt_position(variant.pos))
