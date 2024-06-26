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

from ..uint_range import UIntRange
from ..utils import is_dna, reverse_complement


class DnaStr(str):

    def __init__(self, s: str) -> None:
        if not is_dna(s):
            raise ValueError(f"Invalid DNA sequence: {s}!")
        super().__init__()

    @classmethod
    def parse(cls, s: str | None):
        return cls(s) if s else cls.empty()

    @classmethod
    def empty(cls):
        return cls('')

    def __add__(self, other) -> DnaStr:
        return DnaStr(str(self) + str(other))

    def slice(self, sl: slice) -> DnaStr:
        return DnaStr(self[sl])

    def substr(self, r: UIntRange) -> DnaStr:
        return self.slice(r.to_slice())

    def replace_substr(self, r: UIntRange, alt: str) -> DnaStr:
        assert r.start < len(self) and r.end < len(self)
        return DnaStr(f"{self[:r.start]}{alt}{self[r.end + 1:]}")

    def insert_substr(self, offset: int, alt: str) -> DnaStr:
        """Insert or append a substring"""

        if offset < 0:
            raise ValueError("Invalid string offset: negative!")
        n = len(self)
        if offset > n:
            raise ValueError("Invalid string offset: out of bounds!")
        return DnaStr(
            f"{self[:offset]}{alt}{self[offset:]}" if offset < n else
            self + alt
        )

    def tail(self, n: int) -> DnaStr:
        assert n >= 0 and n <= len(self)
        match n:
            case 0:
                return DnaStr.empty()
            case 1:
                return DnaStr(self[-1])
            case _:
                l: int = len(self)
                return self.substr(UIntRange(l - n - 1 + 1, l))

    def head(self, n: int) -> DnaStr:
        assert n >= 0 and n <= len(self)
        match n:
            case 0:
                return DnaStr.empty()
            case 1:
                return DnaStr(self[0])
            case _:
                return self.substr(UIntRange(0, n - 1))

    @property
    def revc(self):
        return self.__class__(reverse_complement(self))

    def get_at(self, indices: list[int]) -> DnaStr:
        match len(indices):
            case 0:
                return DnaStr.empty()
            case 1:
                return DnaStr(self[indices[0]])
            case _:
                return DnaStr(''.join(self[i] for i in indices))
