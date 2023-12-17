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

from .uint_range import UIntRange


@dataclass(slots=True, frozen=True)
class IntPatternBuilder:
    offset: int
    span: int

    def __post_init__(self) -> None:
        if self.span <= 0:
            raise ValueError(f"Invalid pattern span: {self.span}!")

    def build(self, start: int, length: int) -> list[int]:
        s = start + self.offset
        return list(range(s, s + length + 1 - self.span, self.span))

    def get_range(self, start: int) -> UIntRange:
        return UIntRange.from_length(start, self.span)


pt_codon = IntPatternBuilder(0, 3)
