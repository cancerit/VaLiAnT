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

from __future__ import annotations

from dataclasses import dataclass
from typing import Sized

from .base import PositionRange
from .dna_str import DnaStr


@dataclass(frozen=True)
class Seq(Sized):
    """Non-empty DNA sequence"""

    __slots__ = ['sequence']

    sequence: DnaStr

    def __post_init__(self) -> None:
        if self.sequence.is_empty:
            raise ValueError("Empty DNA sequence!")

    def __len__(self) -> int:
        return len(self.sequence)

    def get_range(self, start: int) -> PositionRange:
        return PositionRange(start, start + len(self) - 1)

    def get_rel_subsequence(self, start: int, pr: PositionRange) -> Seq:
        if pr not in self.get_range(start):
            raise ValueError("Invalid subsequence: out of range!")
        return Seq(self.sequence.slice(pr.start - 1, pr.end))
