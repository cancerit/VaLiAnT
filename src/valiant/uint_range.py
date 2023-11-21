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

from __future__ import annotations

from collections.abc import Container
from dataclasses import dataclass
from typing import Sized, Generic, TypeVar


@dataclass(slots=True, frozen=True)
class UIntRange(Sized, Container):
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError("Invalid range!")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __lt__(self, other) -> bool:
        return (
            self.end < other.end if self.start == other.start else
            self.start < other.start
        )

    def __repr__(self) -> str:
        return f"[{self.start}, {self.end}]"

    def __contains__(self, x) -> bool:
        return self.start <= x <= self.end

    @classmethod
    def from_length(cls, start: int, length: int) -> UIntRange:
        return cls(start, start + length - 1)

    @classmethod
    def from_pos(cls, pos: int) -> UIntRange:
        return UIntRange(pos, pos)

    def to_slice(self, offset: int = 0) -> slice:
        return slice(self.start - offset, self.end - offset + 1)

    def offset(self, offset: int) -> UIntRange:
        return UIntRange(self.start + offset, self.end + offset)


UIntRangeT = TypeVar('UIntRangeT', bound=UIntRange)


@dataclass(slots=True)
class UIntRangeSortedList(Container, Generic[UIntRangeT]):
    ranges: list[UIntRangeT]

    def __post_init__(self) -> None:
        self.ranges.sort()

    def __contains__(self, x) -> bool:
        return any(
            r.start <= x <= r.end
            for r in self.ranges
        )
