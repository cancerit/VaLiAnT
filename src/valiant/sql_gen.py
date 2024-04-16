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
from typing import Any, Iterable

from .db import DbFieldName, DbTableName, cursor


def _sql_boolean_op(op: str, args: list[str], parens: bool = False):
    s = f" {op} "
    x = args if not parens else [f"({x})" for x in args]
    return s.join(x)


def sql_and(args, parens: bool = False) -> str:
    return _sql_boolean_op('and', args, parens=parens)


def sql_or(args, parens: bool = False) -> str:
    return _sql_boolean_op('or', args, parens=parens)


def sql_in_range(
    start_only: bool,
    either: bool = False,
    start_field: DbFieldName = DbFieldName.START,
    end_field: DbFieldName = DbFieldName.END
) -> str:
    assert not (start_only and either)

    return (
        # Either start or end, or both in range
        sql_or([
            sql_and([start_field.ge(), start_field.le()]),
            sql_and([end_field.ge(), end_field.le()])
        ]) if either else
        # Both start and end in range
        sql_and([
            start_field.ge(),
            end_field.le() if not start_only else start_field.le()
        ])
    )


def get_multi_range_check(n: int, start_only: bool = False) -> str:
    return sql_or([sql_in_range(start_only)] * n)


def sql_str(s: str) -> str:
    return f"'{s}'"


def sql_eq_or_in_str_list(field: str, a: list[str]) -> str:
    assert a
    return (
        f"{field} = {sql_str(a[0])}" if len(a) == 1 else
        f"{field} in ({','.join(map(sql_str, a))})"
    )


class SqlQuery(str):

    @classmethod
    def get_delete(cls, t: DbTableName) -> SqlQuery:
        return cls(f"delete from {t.value}")

    @classmethod
    def get_update(cls, t: DbTableName, values: list[DbFieldName], where: str, value: Any = None) -> SqlQuery:
        sets = ','.join(k.sql_eq(value=value) for k in values)
        return cls(f"update {t.value} set {sets} where {where}")

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
    def get_select(cls, t: DbTableName, fields: list[DbFieldName], const: dict[DbFieldName, str] | None = None, where: str | None = None) -> SqlQuery:
        tokens = [
            f.value
            for f in fields
        ] if not const else [
            f.value if f not in const else f"'{const[f]}' as {f.value}"
            for f in fields
        ]

        query = f"select {','.join(tokens)} from {t.value}"

        if where:
            query += f" where {where}"

        return cls(query)

    @classmethod
    def get_select_in_range(
        cls,
        t: DbTableName,
        fields: list[DbFieldName],
        either: bool = False,
        end_field: DbFieldName | None = None,
        start_only: bool = False
    ) -> SqlQuery:
        """
        Select the fields provided where start and end are within a range

        Assumption: the target table has the start and end fields.

        E.g.: select start, ref, alt from variants where start >= ? and end <= ?
        """

        return cls.get_select(t, fields, where=sql_in_range(
            start_only, either=either, end_field=end_field or DbFieldName.END))

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
