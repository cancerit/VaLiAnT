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

from functools import lru_cache

from .dna_str import DnaStr
from ..utils import reverse_complement


class Codon(DnaStr):
    def __init__(self, s: str) -> None:
        if len(s) != 3:
            raise ValueError(f"Invalid codon: {s}!")
        super().__init__(s)

    @lru_cache(maxsize=64)
    def reverse_complement(self) -> Codon:
        return Codon(reverse_complement(self))
