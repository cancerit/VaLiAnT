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

from itertools import chain
from sqlite3 import Connection, Cursor

from .annot_variant import AnnotVariant
from .custom_variant import CustomVariant
from .db import PER_TARGETON_TABLES, cursor, DbTableName, DbFieldName, get_csv_header, select_to_csv
from .exon import Exon
from .experiment_meta import ExperimentMeta
from .oligo_seq import OligoSeq
from .options import Options
from .pam_variant import PamVariant
from .pattern_variant import PatternVariant
from .sql_gen import SqlQuery, SqlScript, get_multi_range_check
from .uint_range import UIntRange
from .utils import bool_to_int_str, get_enum_values
from .variant import RegisteredVariant
from .variant_select import VariantSelectStart


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
            cur.execute(
                sql_insert_ppe_sgrna_ids,
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
from v_background_variants b
where start >= ? and start <= ? and alt_ref_delta != 0
"""




insert_pattern_variant_fields = [
    DbFieldName.POS_A,
    DbFieldName.REF_A,
    DbFieldName.ALT_A,
    DbFieldName.MUTATOR,
    DbFieldName.OLIGO
]


sql_insert_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.PATTERN_VARIANTS, insert_pattern_variant_fields)


def insert_pattern_variants(
    conn: Connection,
    vars: list[OligoSeq[PatternVariant]]
) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_pattern_variants, [
            (
                oligo.variant.pos,
                oligo.variant.ref,
                oligo.variant.alt,
                oligo.variant.mutator,
                oligo.seq
            )
            for oligo in vars
        ])


sql_insert_annot_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.PATTERN_VARIANTS, [
        *insert_pattern_variant_fields,
        DbFieldName.CODON_REF_A,
        DbFieldName.CODON_ALT_A,
        DbFieldName.AA_REF,
        DbFieldName.AA_ALT
    ])


def insert_annot_pattern_variants(conn: Connection, vars: list[OligoSeq[AnnotVariant]]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_annot_pattern_variants, [
            (
                oligo.variant.pos,
                oligo.variant.ref,
                oligo.variant.alt,
                oligo.variant.src,
                oligo.seq,
                oligo.variant.codon_ref,
                oligo.variant.codon_alt,
                oligo.variant.aa_ref,
                oligo.variant.aa_alt
            )
            for oligo in vars
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


clear_per_targeton_tables = SqlScript.from_queries(
    map(SqlQuery.get_delete, PER_TARGETON_TABLES))


select_ppes_in_range = VariantSelectStart.from_table(
    DbTableName.PAM_PROTECTION_EDITS).select_in_range


select_bgs_in_range = VariantSelectStart.from_table(
    DbTableName.BACKGROUND_VARIANTS).select_in_range


select_custom_variants_in_range = VariantSelectStart.from_table(
    DbTableName.CUSTOM_VARIANTS).select_in_range


sql_insert_targeton_ppes = SqlQuery.get_insert_values(
    DbTableName.TARGETON_PAM_PROTECTION_EDITS, [
        DbFieldName.ID,
        DbFieldName.START
    ])


def insert_targeton_ppes(conn: Connection, variants: list[RegisteredVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_targeton_ppes, [
            (v.id, v.pos)
            for v in variants
        ])


def is_table_empty(conn: Connection, t: DbTableName) -> bool:
    query = f"select exists (select 1 from {t.value} limit 1)"
    with cursor(conn) as cur:
        return cur.execute(query).fetchone() == 0


def is_meta_table_empty(conn: Connection) -> bool:
    return is_table_empty(conn, DbTableName.MUTATIONS)


sql_insert_targeton_custom_variants = SqlQuery.get_insert_values(
    DbTableName.TARGETON_CUSTOM_VARIANTS, [
        DbFieldName.ID,
        DbFieldName.START,
        DbFieldName.OLIGO
    ])


def insert_targeton_custom_variants(
    conn: Connection,
    variants: list[OligoSeq[RegisteredVariant]],
    const_regions: list[UIntRange]
) -> None:
    with cursor(conn) as cur:
        # Reference the custom variants in the targeton range
        cur.executemany(sql_insert_targeton_custom_variants, [
            (oligo.variant.id, oligo.variant.pos, oligo.seq)
            for oligo in variants
        ])

        if const_regions:

            # Mark the custom variants that are in constant regions
            query = SqlQuery.get_update(
                DbTableName.TARGETON_CUSTOM_VARIANTS, [DbFieldName.IN_CONST],
                get_multi_range_check(len(const_regions), start_only=True),
                value=1)
            cur.execute(query, tuple(chain.from_iterable(
                r.to_tuple() for r in const_regions)))
