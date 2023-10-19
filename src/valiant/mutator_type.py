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

from enum import Enum
import re

from .utils import parse_opt_int_group


class MutatorType(str, Enum):
    # Parametric span
    DEL = 'del'
    # Single-nucleotide span
    SNV = 'snv'
    # Codon span
    SNV_RE = 'snvre'
    IN_FRAME = 'inframe'
    ALA = 'ala'
    STOP = 'stop'
    AA = 'aa'

    @property
    def is_parametric(self) -> bool:
        return self == MutatorType.DEL


MUTATOR_TYPE_RE = {
    MutatorType.DEL: re.compile(r"^(\d+)del(\d*)$")
}


def parse_mutator_type(s: str) -> tuple[MutatorType, tuple[int, int] | None]:
    for t, pattern in MUTATOR_TYPE_RE.items():
        m = pattern.match(s)
        if m:
            # Parametric pattern
            span = int(m.group(1))
            offset = parse_opt_int_group(m, 2)
            return t, (offset, span)

    # Fixed pattern
    return MutatorType(s), None
