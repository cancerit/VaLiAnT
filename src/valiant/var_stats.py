########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from dataclasses import dataclass, replace
from typing import Iterable

from .uint_range import UIntRange
from .utils import get_end


class OutOfBoundsVar(ValueError):
    def __init__(self, pos: int, r: UIntRange, *args: object) -> None:
        self.pos = pos
        self.r = r
        super().__init__(*args)

    @property
    def msg(self) -> str:
        return f"Out-of-bounds variant: {self.pos} not in {self.r}!"


class OverlappingVar(ValueError):
    pass


@dataclass(slots=True, frozen=True)
class VarStats:
    pos: int
    ref_len: int
    alt_len: int
    alt_pos: int | None = None

    def clone(self):
        return replace(self)

    @property
    def alt_ref_delta(self) -> int:
        return self.alt_len - self.ref_len

    @property
    def ref_end(self) -> int:
        return get_end(self.pos, self.ref_len)

    @property
    def ref_range(self) -> UIntRange:
        return UIntRange.from_length(self.pos, self.ref_len)

    @property
    def alt_range(self) -> UIntRange:
        return UIntRange.from_length(self.pos, self.alt_len)

    def overlaps(self, r: UIntRange) -> bool:
        return self.pos in r or self.ref_end in r


def get_alt_ref_delta(vs: Iterable[VarStats]) -> int:
    return sum(x.alt_ref_delta for x in vs)


def clamp_var_stats_collection(vs: Iterable[VarStats], r: UIntRange) -> list[VarStats]:
    """
    Sort by position and clamp to a range the variant stats

    Fully out-of-bounds variants are not accepted.
    """

    # Safeguard against relative ranges
    assert r.start > 0

    # Sort variant stats by position
    a = sorted(vs, key=lambda x: x.pos)
    if not a:
        return []

    # Test for out of bounds at the start
    v = a[0]
    if v.ref_range not in r:
        raise OutOfBoundsVar(v.pos, r)

    n = len(a)

    if n > 1:

        # Test for overlaps
        if any(a[i].pos == a[i - 1] for i in range(1, n)):
            raise OverlappingVar("Overlapping variants!")

        v = a[-1]

        # Test for out of bounds at the end
        if v.ref_range not in r:
            raise OutOfBoundsVar(v.pos, r)

    return a
