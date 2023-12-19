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

from dataclasses import dataclass, replace
from typing import Sized

from .int_pattern_builder import IntPatternBuilder
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .utils import get_end


@dataclass(slots=True)
class Seq(Sized):
    start: int
    s: DnaStr

    def __len__(self) -> int:
        return len(self.s)

    @property
    def is_length_multiple_of_three(self) -> bool:
        return len(self) % 3 == 0

    @property
    def end(self) -> int:
        return get_end(self.start, len(self))

    def get_range(self) -> UIntRange:
        return UIntRange(self.start, self.end)

    @classmethod
    def from_str(cls, start: int, s: str) -> Seq:
        return cls(start, DnaStr.parse(s))

    def get_rel_range(self, r: UIntRange) -> UIntRange:
        return r.offset(-self.start)

    def _substr(self, s: DnaStr, r: UIntRange, rel: bool = True) -> DnaStr:
        return s.substr(r if rel else self.get_rel_range(r))

    def substr(self, r: UIntRange, rel: bool = True) -> DnaStr:
        # Assumption: if the range is absolute, it is one-based
        # Assumption: if the range is relative, it is zero-based
        return self._substr(self.s, r, rel=rel)

    def replace_substr(self, r: UIntRange, alt: str) -> DnaStr:
        return self.s.replace_substr(self.get_rel_range(r), alt)

    def alter(self, r: UIntRange, alt: str) -> Seq:
        return replace(self, s=self.replace_substr(r, alt))

    def clone(self) -> Seq:
        return replace(self)

    def subseq(self, r: UIntRange, rel: bool = True) -> Seq:
        return Seq(r.start, self.substr(r, rel=rel))

    def subseq_window(self, pt: IntPatternBuilder, start: int = 0, length: int | None = None) -> list[Seq]:
        starts = pt.build(start, (length if length is not None else len(self)) - 1)
        return [
            self.subseq(pt.get_range(start), rel=False)
            for start in starts
        ]

    def subseq_triplets(self) -> list[Seq]:
        pt = IntPatternBuilder(0, 3)
        return self.subseq_window(pt, start=0, length=len(self))

    def get_offset(self, pos: int) -> int:
        return pos - self.start

    def tail(self, n: int) -> DnaStr:
        return self.s.tail(n)

    def head(self, n: int) -> DnaStr:
        return self.s.head(n)
