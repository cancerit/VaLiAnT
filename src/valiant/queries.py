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
from .background_variants import RegisteredBackgroundVariant
from .custom_variant import CustomVariant
from .db import PER_TARGETON_TABLES, cursor, DbTableName, DbFieldName
from .exon import Exon
from .oligo_seq import OligoSeq
from .pam_variant import PamVariant
from .pattern_variant import PatternVariant
from .sql_gen import SqlQuery, SqlScript, get_multi_range_check, sql_eq_or_in_str_list
from .strings.strand import Strand
from .uint_range import UIntRange
from .utils import get_enum_values, safe_group_by
from .variant import RegisteredVariant, Variant
from .variant_group import VariantGroup
from .variant_select import VariantSelectStart, VariantSelectStartEnd


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
        DbFieldName.EXON_INDEX,
        DbFieldName.CDS_PREFIX_LENGTH,
        DbFieldName.CDS_SUFFIX_LENGTH,
        DbFieldName.FIRST_CODON_START
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


def insert_exons(conn: Connection, strand: Strand, exons: list[Exon]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_exon, [
            (
                e.start,
                e.end,
                e.index,
                e.cds_prefix_length,
                e.cds_suffix_length,
                e.get_first_codon_start(strand)
            )
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


# TODO: this shoud NOT assign the exon ID's (those may change per targeton and depend on background!)
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


sql_select_exon_ppes = """
select
    exon_index,
    start,
    ref,
    alt
from v_exon_ppes
"""


def select_exon_ppes(conn: Connection, sgrna_ids: frozenset[str] | None = None) -> dict[int, VariantGroup[Variant]]:
    """Map exon indices to PPE's, optionally filtering by sgRNA ID"""

    query = sql_select_exon_ppes
    if sgrna_ids:
        query = f"{query} where {sql_eq_or_in_str_list('sgrna_id', list(sgrna_ids))}"

    with cursor(conn) as cur:
        res = cur.execute(query).fetchall()
        print(cur.execute("select * from targeton_pam_protection_edits").fetchall())

    return {
        x: VariantGroup.from_variants([
            Variant(*t[1:])
            for t in y
        ])
        for x, y in safe_group_by(res, lambda r: r[0])
    }


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
        DbFieldName.AA_ALT,
        DbFieldName.MUTATION_TYPE
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
                oligo.variant.aa_alt,
                oligo.variant.mutation_type.value
            )
            for oligo in vars
        ])


def insert_enum(conn: Connection, t: DbTableName, cls) -> None:
    with cursor(conn) as cur:
        cur.executemany(SqlQuery.get_insert_names(t), [
            (name,) for name in get_enum_values(cls)
        ])


sql_select_exons_in_range = """
select
    exon_index,
    start,
    end
from exons
where
    (start >= ? and start <= ?) or
    (end >= ? and end <= ?)
"""


def select_exons_in_range(conn: Connection, start: int, end: int) -> list[tuple[int, UIntRange]]:
    with cursor(conn) as cur:
        return [
            (r[0], UIntRange(r[1], r[2]))
            for r in cur.execute(
                sql_select_exons_in_range,
                (start, end, start, end)).fetchall()
        ]


clear_per_targeton_tables = SqlScript.from_queries(
    map(SqlQuery.get_delete, PER_TARGETON_TABLES))


sql_select_ppes_in_range = VariantSelectStart.from_table(
    DbTableName.V_PPE_SGRNA_IDS).query


select_bgs_in_range = VariantSelectStart.from_table(
    DbTableName.BACKGROUND_VARIANTS).select_in_range


select_custom_variants_in_range = VariantSelectStartEnd.from_table(
    DbTableName.V_CUSTOM_VARIANTS).select_in_range


sql_insert_targeton_ppes = SqlQuery.get_insert_values(
    DbTableName.TARGETON_PAM_PROTECTION_EDITS, [
        DbFieldName.ID,
        DbFieldName.START
    ])


# TODO: consider the exons should have background-altered coordinates as well
sql_insert_exon_codon_ppes = """
insert into targeton_exon_codon_ppes (
    ppe_id,
    exon_id,
    codon_index
)
select
    tppe.id,
    e.id, (
        abs(tppe.start - e.first_codon_start) / 3
    ) as codon_index
from targeton_pam_protection_edits tppe
left join pam_protection_edits ppe on ppe.id = tppe.id
left join v_exon_ext e on
    tppe.start >= e.start and
    tppe.start <= e.end;
"""


def insert_targeton_ppes(conn: Connection, variants: list[RegisteredVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_targeton_ppes, [
            (v.id, v.pos)
            for v in variants
        ])

        # Assign exon ID's and codon indices to the PPE's
        cur.execute(sql_insert_exon_codon_ppes)


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


sql_select_ppes_with_offset = """
select
    exon_index,
    ppe_start,
    codon_offset
from v_exon_codon_ppes
where ppe_start >= ? and ppe_start <= ?
order by ppe_start
"""


def select_ppes_with_offset(conn: Connection, r: UIntRange) -> list[tuple[int, int, int]]:
    with cursor(conn) as cur:
        return cur.execute(sql_select_ppes_with_offset, r.to_tuple()).fetchall()


sql_select_background_variants = """
select
    start,
    ref,
    alt,
    id,
    start_exon_index,
    start_codon_index,
    end_exon_index,
    end_codon_index
from v_background_variant_codons
where start >= ? and start <= ?
"""


def select_background_variants(conn: Connection, r: UIntRange) -> list[RegisteredBackgroundVariant]:
    with cursor(conn) as cur:
        return [
            RegisteredBackgroundVariant(*r)
            for r in cur.execute(
                sql_select_background_variants, r.to_tuple()).fetchall()
        ]


sql_select_ppe_bg_codon_overlaps = """
select distinct p.ppe_start
from v_exon_codon_ppes p
left join v_background_variants b on
    p.ppe_start >= b.start and
    p.ppe_start <= b.ref_end
where
    p.exon_index is not null and
    b.id is not null
"""


def select_ppe_bg_codon_overlaps(conn: Connection) -> list[int]:
    with cursor(conn) as cur:
        return [
            r[0]
            for r in cur.execute(sql_select_ppe_bg_codon_overlaps).fetchall()
        ]
