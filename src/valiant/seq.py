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

from dataclasses import dataclass
from typing import Sized

from .strings.dna_str import DnaStr
from .uint_range import UIntRange


@dataclass(slots=True)
class Seq(Sized):
    start: int
    s: DnaStr

    def __len__(self) -> int:
        return len(self.s)

    @property
    def end(self) -> int:
        return self.start + max(0, len(self) - 1)

    def get_range(self) -> UIntRange:
        return UIntRange(self.start, self.end)

    @classmethod
    def from_str(cls, start: int, s: str) -> Seq:
        return cls(start, DnaStr.parse(s))

    def substr(self, r: UIntRange, rel: bool = True) -> DnaStr:
        return self.s.substr(
            r if rel else
            r.offset(-self.start)
        )

    def subseq(self, r: UIntRange, rel: bool = True) -> Seq:
        return Seq(r.start, self.substr(r, rel=rel))
