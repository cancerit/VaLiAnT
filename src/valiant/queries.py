########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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
from sqlite3 import Connection, Cursor, IntegrityError

from .annot_variant import AnnotVariant
from .custom_variant import CustomVariant
from .db import PER_TARGETON_TABLES, VARIANT_FIELDS, cursor, DbTableName, DbFieldName, PER_CONTIG_TABLES
from .exon import Exon
from .oligo_seq import OligoSeq
from .pam_variant import InvalidPamVariant, PamVariant
from .sql_gen import SqlQuery, SqlScript, get_multi_range_check, sql_eq_or_in_str_list
from .strings.dna_str import DnaStr
from .strings.strand import Strand
from .uint_range import UIntRange
from .variant import PatternVariant, RegisteredVariant, VarStats, Variant
from .variant_select import VariantSelectStartEnd


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


insert_pattern_variant_fields = [
    DbFieldName.POS_R,
    DbFieldName.POS_A,
    DbFieldName.REF_A,
    DbFieldName.ALT_A,
    DbFieldName.MUTATOR,
    DbFieldName.OLIGO
]


sql_insert_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.ALT_PATTERN_VARIANTS, insert_pattern_variant_fields)


def insert_pattern_variants(
    conn: Connection,
    vars: list[OligoSeq[PatternVariant]]
) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_pattern_variants, [
            (
                oligo.ref_start,
                oligo.variant.pos,
                oligo.variant.ref,
                oligo.variant.alt,
                oligo.variant.mutator,
                oligo.seq
            )
            for oligo in vars
        ])


sql_insert_annot_pattern_variants = SqlQuery.get_insert_values(
    DbTableName.ALT_PATTERN_VARIANTS, [
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
                oligo.ref_start,
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


sql_select_exons_in_range = """
select
    exon_index,
    start,
    end
from exons
where
    (? >= start and ? <= end) or
    (? >= start and ? <= end)
"""


def select_exons_in_range(conn: Connection, start: int, end: int) -> list[tuple[int, UIntRange]]:
    with cursor(conn) as cur:
        return [
            (r[0], UIntRange(r[1], r[2]))
            for r in cur.execute(
                sql_select_exons_in_range,
                (start, start, end, end)).fetchall()
        ]


clear_per_contig_tables = SqlScript.from_queries(
    map(SqlQuery.get_delete, PER_CONTIG_TABLES))


clear_per_targeton_tables = SqlScript.from_queries(
    map(SqlQuery.get_delete, PER_TARGETON_TABLES))


select_custom_variants_in_range = VariantSelectStartEnd.from_table(
    DbTableName.V_CUSTOM_VARIANTS).select_in_range


sql_insert_targeton_ppes = SqlQuery.get_insert_values(
    DbTableName.TARGETON_PAM_PROTECTION_EDITS, [
        DbFieldName.ID,
        DbFieldName.START
    ])


# Consider the exons should have background-altered coordinates as well
# TODO: evaluate implications of discarding PPE's in noncoding regions vs.
#  making the exon ID and codon index nullable
#  (see also the view that references this table)
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
    tppe.start <= e.end
where e.id is not null;
"""


def insert_targeton_ppes(conn: Connection, variants: list[RegisteredVariant]) -> None:
    with cursor(conn) as cur:
        cur.executemany(sql_insert_targeton_ppes, [
            (v.id, v.pos)
            for v in variants
        ])

        try:

            # Assign exon ID's and codon indices to the PPE's
            cur.execute(sql_insert_exon_codon_ppes)

        except IntegrityError as ex:
            if ex.sqlite_errorname == 'SQLITE_CONSTRAINT_PRIMARYKEY':
                raise InvalidPamVariant


def is_table_empty(conn: Connection, t: DbTableName) -> bool:
    query = f"select exists (select 1 from {t.value} limit 1)"
    with cursor(conn) as cur:
        return cur.execute(query).fetchone()[0] == 0


def is_meta_table_empty(conn: Connection) -> bool:
    return is_table_empty(conn, DbTableName.V_META)


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
    ppe_ref_start,
    ppe_start,
    codon_offset
from v_exon_codon_ppes
where ppe_start >= ? and ppe_start <= ?
order by ppe_start
"""


def select_ppes_with_offset(conn: Connection, r: UIntRange) -> list[tuple[int, int, int, int]]:
    with cursor(conn) as cur:
        return cur.execute(sql_select_ppes_with_offset, r.to_tuple()).fetchall()


sql_select_background_variants = """
select
    start,
    ref,
    alt,
    id
from background_variants
where start >= ? and start <= ?
order by start
"""


def select_background_variants(conn: Connection, r: UIntRange) -> list[RegisteredVariant]:
    with cursor(conn) as cur:
        return [
            RegisteredVariant(*r)
            for r in cur.execute(
                sql_select_background_variants, r.to_tuple()).fetchall()
        ]


sql_select_overlapping_background_variants = """
select start, ref, alt
from v_background_variants
where
    (start >= ? and start <= ?) or
    (ref_end >= ? and ref_end <= ?)
order by start, ref, alt
"""


def select_overlapping_background_variants(conn: Connection, r: UIntRange) -> list[Variant]:
    t = r.to_tuple()
    with cursor(conn) as cur:
        return [
            Variant(*r)
            for r in cur.execute(
                sql_select_overlapping_background_variants, (*t, *t)).fetchall()
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


sql_select_background_variant_stats = SqlQuery.get_select_in_range(
    DbTableName.V_BACKGROUND_VARIANTS, [
        DbFieldName.START,
        DbFieldName.REF_LEN,
        DbFieldName.ALT_LEN
    ],
    either=True,
    end_field=DbFieldName.REF_END)


def select_background_variant_stats(conn: Connection, r: UIntRange) -> list[VarStats]:
    with cursor(conn) as cur:
        return [
            VarStats(*row)
            for row in cur.execute(
                sql_select_background_variant_stats,
                (r.start, r.end, r.start, r.end))
        ]


def _get_select_ppes_by_sgrna_id(sgrna_ids: frozenset[str]) -> SqlQuery:
    """
    Generate a query to select PPE's by their sgRNA ID's

    E.g. (formatting aside):
        select start, ref, alt, id
        from v_ppe_sgrna_ids
        where sgrna_id in ('sgRNA_a','sgRNA_b')
    """

    return SqlQuery.get_select(
        DbTableName.V_PPE_SGRNA_IDS,
        VARIANT_FIELDS,
        where=sql_eq_or_in_str_list(
            'sgrna_id', list(sgrna_ids)))


def select_ppes_by_sgrna_id(conn: Connection, sgrna_ids: frozenset[str]) -> list[RegisteredVariant]:
    query = _get_select_ppes_by_sgrna_id(sgrna_ids)
    with cursor(conn) as cur:
        return [
            RegisteredVariant(r[0], DnaStr(r[1]), DnaStr(r[2]), r[3])
            for r in cur.execute(query)
        ]
