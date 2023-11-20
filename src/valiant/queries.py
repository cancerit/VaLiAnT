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

from sqlite3 import Connection, Cursor

from .custom_variant import CustomVariant
from .db import cursor, get_sql_insert, get_sql_insert_values, DbTableName, DbFieldName, get_sql_select_name
from .exon import Exon
from .mutator_type import MutatorType
from .pam_variant import PamVariant
from .pattern_variant import PatternVariant
from .utils import get_enum_values


class NoRowId(Exception):
    pass


class SqlInsertFailed(Exception):
    pass


def get_sql_insert_names(t: DbTableName) -> str:
    return get_sql_insert_values(t, [DbFieldName.NAME])


sql_insert_custom_variant_collection = get_sql_insert_names(
    DbTableName.CUSTOM_VARIANT_COLLECTIONS)

sql_insert_custom_variant: str = get_sql_insert_values(
    DbTableName.CUSTOM_VARIANTS,
    [
        DbFieldName.COLLECTION_ID,
        DbFieldName.VAR_ID,
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT,
        DbFieldName.VCF_NT
    ])

sql_insert_exon = get_sql_insert_values(
    DbTableName.EXONS,
    [
        DbFieldName.START,
        DbFieldName.END,
        DbFieldName.EXON_INDEX
    ])

sql_insert_background_variant = get_sql_insert_values(
    DbTableName.BACKGROUND_VARIANTS,
    [
        DbFieldName.VAR_ID,
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT
    ])

sql_insert_ppe = get_sql_insert_values(
    DbTableName.PAM_PROTECTION_EDITS,
    [
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT
    ])

sql_insert_ppe_sgrna_ids = get_sql_insert_values(
    DbTableName.PAM_PROTECTION_EDIT_SGRNA_IDS,
    [
        DbFieldName.VAR_PPE_ID,
        DbFieldName.SGRNA_ID
    ])

sql_insert_sgrna_id = get_sql_insert_names(DbTableName.SGRNA_IDS)
sql_insert_mutator_types = get_sql_insert_names(DbTableName.MUTATOR_TYPES)


def insert_custom_variant_collection(conn: Connection, name: str, vars: list[CustomVariant]) -> None:
    with cursor(conn) as cur:

        # Insert the collection
        cur.execute(sql_insert_custom_variant_collection, (name,))
        collection_id = cur.lastrowid
        assert collection_id is not None

        # Insert the variants
        cur.executemany(sql_insert_custom_variant, [
            (collection_id, v.id, v.pos, v.ref, v.alt, v.vcf_nt or '')
            for v in vars
        ])


def insert_exons(conn: Connection, exons: list[Exon]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_exon, [
            (e.start, e.end, e.index)
            for e in exons
        ])


def insert_background_variants(conn: Connection, vars: list[CustomVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_background_variant, [
            (v.id, v.pos, v.ref, v.alt)
            for v in vars
        ])


def _insert_returning_row_id(cur: Cursor, query: str, params: tuple) -> int:
    assert query.lower().startswith('insert into ')
    cur.execute(query, params)
    row_id = cur.lastrowid
    if row_id is None:
        raise NoRowId
    return row_id


def insert_pam_protection_edits(conn: Connection, vars: list[PamVariant]) -> None:
    def get_sgrna_id(c: Cursor, name: str) -> int:

        # Insert sgRNA ID
        try:
            row_id = _insert_returning_row_id(c, sql_insert_sgrna_id, (name,))
        except NoRowId:
            raise SqlInsertFailed("Failed to insert PPE!")
        return row_id

    sgrna_ids = {v.sgrna_id for v in vars}
    with cursor(conn) as cur:

        # Insert all sgRNA ID's
        # Assumption: this is run only once
        #  (otherwise, duplicates would need handling)
        sgrna_name_ids: dict[str, int] = {
            name: get_sgrna_id(cur, name)
            for name in sorted(sgrna_ids)
        }

        for v in vars:

            # Insert variant
            cur.execute(sql_insert_ppe, (v.pos, v.ref, v.alt))
            var_id = cur.lastrowid
            assert var_id is not None

            # Match variant and sgRNA ID
            cur.execute(sql_insert_ppe_sgrna_ids,
                        (var_id, sgrna_name_ids[v.sgrna_id]))


# TODO: consider start positions get offset as they go...
sql_select_offsets = """
select
    start,
    sum(alt_ref_delta) over (
        order by start
        rows between unbounded preceding and current row
    ) as offset
from background_variants_v b
where start >= ? and start <= ? and alt_ref_delta != 0
"""

sql_insert_offsets = get_sql_insert(
    DbTableName.BACKGROUND_OFFSETS,
    [DbFieldName.REF_POS, DbFieldName.OFFSET],
    sql_select_offsets)


def insert_gene_offsets(conn: Connection, start: int, end: int) -> None:
    with cursor(conn) as cur:
        cur.execute(sql_insert_offsets, (start, end))


sql_insert_pattern_variants = get_sql_insert_values(
    DbTableName.PATTERN_VARIANTS, [
        DbFieldName.POS_A,
        DbFieldName.REF_A,
        DbFieldName.ALT_A,
        DbFieldName.MUTATOR
    ])


sql_select_mutator_type_id = get_sql_select_name(DbTableName.MUTATOR_TYPES)


def insert_pattern_variants(conn: Connection, vars: list[PatternVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_pattern_variants, [
            (v.pos, v.ref, v.alt, v.mutator)
            for v in vars
        ])


def insert_enum(conn: Connection, t: DbTableName, cls) -> None:
    with cursor(conn) as cur:
        cur.executemany(get_sql_insert_names(t), [
            (name,) for name in get_enum_values(cls)
        ])


def insert_mutator_types(conn: Connection) -> None:
    insert_enum(conn, DbTableName.MUTATOR_TYPES, MutatorType)


sql_select_exons_in_range = """
select exon_index
from exons
where
    (start >= ? and start <= ?) or
    (end >= ? and end <= ?)
"""


def select_exons_in_range(conn: Connection, start: int, end: int) -> list[int]:
    with cursor(conn) as cur:
        return [
            r[0]
            for r in cur.execute(
                sql_select_exons_in_range,
                (start, end, start, end)).fetchall()
        ]
