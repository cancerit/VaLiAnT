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

import sqlite3
from contextlib import contextmanager
from enum import Enum
from sqlite3 import Connection, Cursor
from typing import Generator

from .enums import DbTable


class DbTableName(str, Enum):
    PAM_PROTECTION_EDITS = 'pam_protection_edits'
    PAM_PROTECTION_EDIT_SGRNA_IDS = 'pam_protection_edit_sgrna_ids'
    BACKGROUND_VARIANTS = 'background_variants'
    EXONS = 'exons'
    CUSTOM_VARIANTS = 'custom_variants'
    CUSTOM_VARIANT_COLLECTIONS = 'custom_variant_collections'
    SGRNA_IDS = 'sgrna_ids'
    BACKGROUND_OFFSETS = 'background_offsets'
    PATTERN_VARIANTS = 'pattern_variants'


class DbFieldName(str, Enum):
    START = 'start'
    END = 'end'
    NAME = 'name'
    VAR_PPE_ID = 'var_ppe_id'
    SGRNA_ID = 'sgrna_id'
    COLLECTION_ID = 'collection_id'
    EXON_INDEX = 'exon_index'
    REF = 'ref'
    ALT = 'alt'
    ID = 'id'
    VAR_ID = 'var_id'
    VCF_NT = 'vcf_nt'
    REF_POS = 'ref_pos'
    OFFSET = 'offset'
    POS_R = 'pos_r'
    REF_R = 'ref_r'
    ALT_R = 'alt_r'
    POS_A = 'pos_a'
    REF_A = 'ref_a'
    ALT_A = 'alt_a'
    MUTATOR = 'mutator'
    CODON_REF_A = 'codon_ref_a'
    CODON_ALT_A = 'codon_alt_a'
    AA_REF = 'aa_ref'
    AA_ALT = 'aa_alt'
    # MUTATOR_TYPE_ID = 'mutator_type_id'


@contextmanager
def get_db_conn() -> Generator[sqlite3.Connection, None, None]:
    with sqlite3.connect('valiant.db' or ':memory:') as conn:
        try:
            cur = conn.cursor()
            yield conn
        finally:
            pass
            # conn.close()


def get_insert(table_name: DbTable, fields: list[str], values: list[str]) -> str:
    fields_s = ','.join(fields)
    values_s = ','.join(f"({v})" for v in values)
    return f"insert into {table_name.value}({fields_s})values{values_s}"


def get_create_view(name: str, fields: list[str], select: str) -> str:
    fields_s = ','.join(fields)
    return f"create view {name}({fields_s}) as {select}"


def create_range_index(table: DbTable) -> str:
    table_name = table.name
    start = DbFieldName.START.value
    end = DbFieldName.END.value
    index_name: str = f"{table_name}_{start}_{end}_idx"
    return f"create index {index_name} on {table_name}({start}, {end})"


def get_insert_ranges(table_name: DbTable, ranges: list[tuple[int, int]]) -> str:
    start = DbFieldName.START.value
    end = DbFieldName.END.value
    return get_insert(table_name, [start, end], [f"{x},{y}" for x, y in ranges])


@contextmanager
def cursor(conn: Connection) -> Generator[Cursor, None, None]:
    cur = conn.cursor()
    try:
        yield cur
    finally:
        cur.close()


def get_sql_insert(t: DbTableName, fields: list[DbFieldName], values: str) -> str:
    return f"insert into {t.value}({','.join(f.value for f in fields)}){values}"


def get_sql_select_name(t: DbTableName) -> str:
    return f"select id from {t.value} where name = ? limit 1"


def get_sql_insert_values(t: DbTableName, fields: list[DbFieldName], fks: dict[DbFieldName, DbTableName] | None = None) -> str:
    """
    Generate a SQL insert statement

    E.g.: insert into t (x, y) values (?, ?)
    """

    values = [
        f"({get_sql_select_name(fks[f])})" if f in fks else '?'
        for f in fields
    ] if fks else '?' * len(fields)

    return get_sql_insert(t, fields, f"values({','.join(values)})")


def dump_table(conn: Connection, t: DbTableName, fh) -> None:
    for r in conn.execute(f"select * from {t.value}"):
        fh.write('\t'.join(str(x) if x is not None else '.' for x in r))
        fh.write('\n')


def dump_all(conn: Connection) -> None:
    for t in DbTableName:
        with open(f"{t.value}_dump.tsv", 'w') as fh:
            dump_table(conn, t, fh)
