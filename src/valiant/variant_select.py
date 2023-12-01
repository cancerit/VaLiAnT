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

from dataclasses import dataclass
from sqlite3 import Connection

from .db import VARIANT_FIELDS, DbTableName, cursor
from .sql_gen import SqlQuery
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .variant import Variant


@dataclass(slots=True)
class VariantSelect:
    query: SqlQuery

    @classmethod
    def from_table(cls, t: DbTableName) -> VariantSelect:
        # TODO: consider partially overlapping variant filtering
        return cls(SqlQuery.get_select_in_range(t, VARIANT_FIELDS))

    def select_in_range(self, conn: Connection, r: UIntRange) -> list[Variant]:
        with cursor(conn) as cur:
            return [
                Variant(r[0], DnaStr(r[1]), DnaStr(r[2]))
                for r in cur.execute(self.query, (r.start, r.end)).fetchall()
            ]
