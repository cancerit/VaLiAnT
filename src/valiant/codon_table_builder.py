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

from .codon_table_loader import load_codon_table_rows
from .codon_table_row import CodonTableRow
from .codon_table import CodonTable
from .strings.strand import Strand


@dataclass(slots=True, frozen=True)
class CodonTableBuilder:
    rows: list[CodonTableRow]

    @classmethod
    def load(cls, fp: str):
        return cls(load_codon_table_rows(fp))

    def build(self, strand: Strand) -> CodonTable:
        return CodonTable.from_list(
            self.rows, rc=strand.is_minus)
