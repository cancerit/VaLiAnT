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

from .annot_variant import AnnotVariant
from .custom_variant import CustomVariant
from .db import PER_TARGETON_TABLES, cursor, DbTableName, DbFieldName, get_csv_header, select_to_csv
from .experiment_meta import ExperimentMeta
from .exon import Exon
from .pam_variant import PamVariant
from .pattern_variant import PatternVariant
from .sql_gen import SqlQuery, SqlScript
from .utils import get_enum_values


class NoRowId(Exception):
    pass


class SqlInsertFailed(Exception):
    pass


sql_insert_custom_variant_collection = SqlQuery.get_insert_names(
    DbTableName.CUSTOM_VARIANT_COLLECTIONS)

sql_insert_custom_variant: str = SqlQuery.get_insert_values(
    DbTableName.CUSTOM_VARIANTS,
    [
        DbFieldName.COLLECTION_ID,
        DbFieldName.VAR_ID,
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT,
        DbFieldName.VCF_NT
    ])

sql_insert_exon = SqlQuery.get_insert_values(
    DbTableName.EXONS,
    [
        DbFieldName.START,
        DbFieldName.END,
        DbFieldName.EXON_INDEX
    ])

sql_insert_background_variant = SqlQuery.get_insert_values(
    DbTableName.BACKGROUND_VARIANTS,
    [
        DbFieldName.VAR_ID,
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT
    ])

sql_insert_ppe = SqlQuery.get_insert_values(
    DbTableName.PAM_PROTECTION_EDITS,
    [
        DbFieldName.START,
        DbFieldName.REF,
        DbFieldName.ALT
    ])

sql_insert_ppe_sgrna_ids = SqlQuery.get_insert_values(
    DbTableName.PAM_PROTECTION_EDIT_SGRNA_IDS,
    [
        DbFieldName.VAR_PPE_ID,
        DbFieldName.SGRNA_ID
    ])

sql_insert_sgrna_id = SqlQuery.get_insert_names(DbTableName.SGRNA_IDS)


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


sql_insert_exon_codon_ppes = """
insert into exon_codon_ppes (
    ppe_id,
    exon_id,
    codon_index
)
select
    p.id,
    e.id, (
        (p.start - e.start - e.cds_prefix_length) / 3
    ) as codon_index
from pam_protection_edits p
left join v_exon_ext e on
    p.start >= e.start and
    p.start <= e.end;
"""


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

        # Assign exon ID's and codon indices to the PPE's
        cur.execute(sql_insert_exon_codon_ppes)


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





insert_pattern_variant_fields = [
    DbFieldName.POS_A,
    DbFieldName.REF_A,
    DbFieldName.ALT_A,
    DbFieldName.MUTATOR
]


sql_insert_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.PATTERN_VARIANTS, insert_pattern_variant_fields)


def insert_pattern_variants(
    conn: Connection,
    vars: list[PatternVariant]
) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_pattern_variants, [
            (v.pos, v.ref, v.alt, v.mutator)
            for v in vars
        ])


sql_insert_annot_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.PATTERN_VARIANTS, [
        *insert_pattern_variant_fields,
        DbFieldName.CODON_REF_A,
        DbFieldName.CODON_ALT_A,
        DbFieldName.AA_REF,
        DbFieldName.AA_ALT
    ])


def insert_annot_pattern_variants(conn: Connection, vars: list[AnnotVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_annot_pattern_variants, [
            (
                v.pos, v.ref, v.alt, v.src,
                v.codon_ref, v.codon_alt,
                v.aa_ref, v.aa_alt
            )
            for v in vars
        ])


def insert_enum(conn: Connection, t: DbTableName, cls) -> None:
    with cursor(conn) as cur:
        cur.executemany(SqlQuery.get_insert_names(t), [
            (name,) for name in get_enum_values(cls)
        ])


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


select_meta_fields = [
    DbFieldName.SPECIES,
    DbFieldName.ASSEMBLY,
    DbFieldName.REF_START,
    DbFieldName.REF,
    DbFieldName.ALT,
    DbFieldName.REF_AA,
    DbFieldName.ALT_AA,
    DbFieldName.VCF_VAR_ID,
    DbFieldName.VCF_ALIAS,
    DbFieldName.MUTATOR,
    DbFieldName.SGRNA_IDS
]
select_meta_header = get_csv_header(select_meta_fields)


def dump_metadata(conn: Connection, exp: ExperimentMeta, fp: str) -> None:
    sql_select_meta = SqlQuery.get_select(
        DbTableName.V_META,
        select_meta_fields,
        const={
            DbFieldName.SPECIES: exp.species,
            DbFieldName.ASSEMBLY: exp.assembly
        })

    select_to_csv(
        conn,
        sql_select_meta,
        len(select_meta_fields),
        select_meta_header,
        fp)


clear_per_targeton_tables = SqlScript.from_queries(
    map(SqlQuery.get_delete, PER_TARGETON_TABLES))
