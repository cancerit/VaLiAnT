########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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
from enum import IntEnum

from .csv import load_csv
from ..uint_range import UIntRange


class BedField(IntEnum):
    CONTIG = 0
    START = 1
    END = 2


@dataclass(slots=True)
class BedLoader:
    contig: str

    def load(self, fp: str) -> list[UIntRange]:
        return [
            UIntRange(
                int(r[BedField.START]) + 1,
                int(r[BedField.END]))
            for r in load_csv(fp, delimiter='\t')
            if r[BedField.CONTIG] == self.contig
        ]
