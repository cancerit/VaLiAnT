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

from ..constants import STOP
from ..enums import MutationType


class TranslationSymbol(str):
    def __init__(self, s: str) -> None:
        if len(s) != 1 and s != STOP:
            raise ValueError(f"Invalid translation symbol: {s}!")
        super().__init__()

    @lru_cache(maxsize=512)
    def get_aa_change(self, aa_b: TranslationSymbol) -> MutationType:
        return (
            MutationType.NONSENSE if aa_b == STOP else
            MutationType.SYNONYMOUS if aa_b == self else
            MutationType.MISSENSE
        )
