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

from dataclasses import dataclass
from collections.abc import Container, Sized


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
            raise ValueError("Invalid range!")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __contains__(self, b) -> bool:
        if isinstance(b, int):
            return b >= self.start and b <= self.end
        raise TypeError("Unsupported operand type!")

    def __eq__(self, b) -> bool:
        if isinstance(b, UIntRange):
            return self.start == b.start and self.end == b.end
        raise TypeError("Unsupported operand type!")

    def __sub__(self, b) -> 'UIntRange':
        if isinstance(b, int):
            return UIntRange(self.start - b, self.end - b)
        raise TypeError("Unsupported operand type!")

    def to_slice(self) -> slice:
        return slice(self.start, self.end + 1)
