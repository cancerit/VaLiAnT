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

from collections.abc import Container, Sized
from dataclasses import dataclass
from typing import List, NoReturn, TypeVar


UIntRangeT = TypeVar('UIntRangeT', bound='UIntRange')


def err_unsupported_operand() -> NoReturn:
    raise TypeError("Unsupported operand type!")


@dataclass(frozen=True)
class UIntRange(Sized, Container):
    """End-inclusive integer range"""

    __slots__ = ['start', 'end']

    start: int
    end: int

    def __post_init__(self) -> None:
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise TypeError("Invalid range boundary types!")
        if self.start < 0:
            raise ValueError("Invalid range start!")
        if self.end < 0:
            raise ValueError("Invalid range end!")
        if self.end < self.start:
            raise ValueError(f"Invalid range [{self.start}, {self.end}]!")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __contains__(self, b) -> bool:
        if isinstance(b, int):
            return self.start <= b <= self.end
        elif isinstance(b, UIntRange):
            return b.start >= self.start and b.end <= self.end
        err_unsupported_operand()

    def __eq__(self, b) -> bool:
        if isinstance(b, UIntRange):
            return self.start == b.start and self.end == b.end
        err_unsupported_operand()

    def __lt__(self, other) -> bool:
        return (
            self.end < other.end if self.start == other.start else
            self.start < other.start
        )

    def __add__(self, b) -> UIntRange:
        if isinstance(b, int):
            return UIntRange(self.start + b, self.end + b)
        err_unsupported_operand()

    def __sub__(self, b) -> UIntRange:
        if isinstance(b, int):
            return UIntRange(self.start - b, self.end - b)
        err_unsupported_operand()

    def to_slice(self) -> slice:
        return slice(self.start, self.end + 1)

    def _from_uintr(self, x: UIntRange) -> UIntRange:
        return x

    def diff(self, others: List[UIntRangeT], skip_out_of_range: bool = False) -> List[UIntRangeT]:
        """Generate the difference set of ranges"""

        new_start: int = self.start
        result = []

        for other in sorted(others):
            if other not in self:
                if skip_out_of_range:
                    continue
                else:
                    raise ValueError("Range out of range!")
            if other.start > new_start:
                result.append(self._from_uintr(UIntRange(new_start, other.start - 1)))
            new_start = other.end + 1

        if new_start <= self.end:
            result.append(self._from_uintr(UIntRange(new_start, self.end)))

        return result
