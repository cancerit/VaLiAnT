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

from .strings.codon import Codon
from .strings.translation_symbol import TranslationSymbol
from .utils import reverse_complement


@dataclass(slots=True, frozen=True)
class CodonTableRow:
    codon: Codon
    aa: TranslationSymbol
    rank: int

    def reverse_complement(self) -> CodonTableRow:
        return replace(self, codon=reverse_complement(self.codon))
