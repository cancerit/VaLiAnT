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

from sqlite3 import Connection
from typing import Iterable

from .db import DbFieldName, DbTableName, cursor


class SqlQuery(str):

    @classmethod
    def get_delete(cls, t: DbTableName) -> SqlQuery:
        return cls(f"delete from {t.value}")

    @classmethod
    def get_insert(cls, t: DbTableName, fields: list[DbFieldName], values: str) -> SqlQuery:
        return cls(f"insert into {t.value}({','.join(f.value for f in fields)}){values}")

    @classmethod
    def get_insert_values(cls, t: DbTableName, fields: list[DbFieldName], fks: dict[DbFieldName, DbTableName] | None = None) -> SqlQuery:
        """
        Generate a SQL insert statement

        E.g.: insert into t (x, y) values (?, ?)
        """

        values = [
            f"({cls.get_select_name(fks[f])})" if f in fks else '?'
            for f in fields
        ] if fks else '?' * len(fields)

        return cls.get_insert(t, fields, f"values({','.join(values)})")

    @classmethod
    def get_insert_names(cls, t: DbTableName) -> SqlQuery:
        return cls.get_insert_values(t, [DbFieldName.NAME])

    @classmethod
    def get_select(cls, t: DbTableName, fields: list[DbFieldName], const: dict[DbFieldName, str] | None = None) -> SqlQuery:
        tokens = [
            f.value
            for f in fields
        ] if not const else [
            f.value if f not in const else f"'{const[f]}' as {f.value}"
            for f in fields
        ]

        return cls(f"select {','.join(tokens)} from {t.value}")

    @classmethod
    def get_select_name(cls, t: DbTableName) -> SqlQuery:
        return cls(f"select id from {t.value} where name = ? limit 1")


class SqlScript(SqlQuery):

    @classmethod
    def from_queries(cls, queries: Iterable[SqlQuery]) -> SqlScript:
        return cls(';'.join(queries) + ';')

    def execute(self, conn: Connection) -> None:
        with cursor(conn) as cur:
            cur.executescript(self)
