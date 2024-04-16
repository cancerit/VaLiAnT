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

from abc import ABC
from dataclasses import dataclass
from sqlite3 import Connection
from typing import ClassVar

from .db import VARIANT_FIELDS, DbTableName, cursor
from .sql_gen import SqlQuery
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .variant import RegisteredVariant


@dataclass(slots=True)
class VariantSelect(ABC):
    _START_ONLY: ClassVar[bool]

    query: SqlQuery

    @classmethod
    def from_table(cls, t: DbTableName):
        return cls(SqlQuery.get_select_in_range(
            t, VARIANT_FIELDS, start_only=cls._START_ONLY))

    def select_in_range(self, conn: Connection, r: UIntRange) -> list[RegisteredVariant]:
        with cursor(conn) as cur:
            return [
                RegisteredVariant(r[0], DnaStr(r[1]), DnaStr(r[2]), r[3])
                for r in cur.execute(self.query, r.to_tuple()).fetchall()
            ]


@dataclass(slots=True)
class VariantSelectStartEnd(VariantSelect):
    _START_ONLY = False
    # TODO: consider partially overlapping variant filtering
