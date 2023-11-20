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

from .variant import Variant


@dataclass
class PatternVariant(Variant):
    mutator: str

    @classmethod
    def from_variant(cls, mutator: str, v: Variant) -> PatternVariant:
        return cls(pos=v.pos, ref=v.ref, alt=v.alt, mutator=mutator)

    def ltrim(self, n: int) -> PatternVariant:
        """
        Trim the first n nucleotides from both REF and ALT and offset
        the start position by the same amount

        Only valid for codon replacements.
        """

        return replace(
            self,
            pos=self.pos + n,
            ref=self.ref.ltrim(n),
            alt=self.alt.ltrim(n))