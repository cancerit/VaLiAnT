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

from ..utils import is_dna


class DnaStr(str):
    """Non-ambiguous DNA string"""

    def __init__(self, s: str) -> None:
        super().__init__()
        if self and not is_dna(self):
            raise ValueError(f"Invalid DNA sequence '{self}'!")

    @property
    def is_empty(self) -> bool:
        return len(self) == 0

    @classmethod
    def empty(cls) -> DnaStr:
        return cls('')

    def slice(self, start: int, end: int) -> DnaStr:
        return DnaStr(self[start:end])
