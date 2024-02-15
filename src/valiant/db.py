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
from typing import Any, Generator

from .utils import get_ddl_path


class DbTableName(str, Enum):
    PAM_PROTECTION_EDITS = 'pam_protection_edits'
    TARGETON_PAM_PROTECTION_EDITS = 'targeton_pam_protection_edits'
    PAM_PROTECTION_EDIT_SGRNA_IDS = 'pam_protection_edit_sgrna_ids'
    BACKGROUND_VARIANTS = 'background_variants'
    V_BACKGROUND_VARIANTS = 'v_background_variants'
    EXONS = 'exons'
    CUSTOM_VARIANTS = 'custom_variants'
    CUSTOM_VARIANT_COLLECTIONS = 'custom_variant_collections'
    SGRNA_IDS = 'sgrna_ids'
    ALT_PATTERN_VARIANTS = 'alt_pattern_variants'
    V_META = 'v_meta'
    MUTATIONS = 'mutations'
    TARGETON_CUSTOM_VARIANTS = 'targeton_custom_variants'
    V_CUSTOM_VARIANTS = 'v_custom_variants'
    V_PPE_SGRNA_IDS = 'v_ppe_sgrna_ids'
    TARGETON_EXON_CODON_PPES = 'targeton_exon_codon_ppes'


# TODO: Verify the custom variant collections can be loaded independently...
PER_CONTIG_TABLES: list[DbTableName] = [
    DbTableName.PAM_PROTECTION_EDIT_SGRNA_IDS,
    DbTableName.PAM_PROTECTION_EDITS,
    DbTableName.BACKGROUND_VARIANTS,
    DbTableName.CUSTOM_VARIANT_COLLECTIONS,
    DbTableName.CUSTOM_VARIANTS
]


PER_TARGETON_TABLES: set[DbTableName] = {
    DbTableName.MUTATIONS,
    DbTableName.ALT_PATTERN_VARIANTS,
    DbTableName.TARGETON_PAM_PROTECTION_EDITS,
    DbTableName.TARGETON_CUSTOM_VARIANTS,
    DbTableName.TARGETON_EXON_CODON_PPES
}


class DbFieldName(str, Enum):
    START = 'start'
    END = 'end'
    REF_END = 'ref_end'
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
    REF_START = 'ref_start'
    REF_AA = 'ref_aa'
    ALT_AA = 'alt_aa'
    VCF_VAR_ID = 'vcf_var_id'
    VCF_ALIAS = 'vcf_alias'
    SGRNA_IDS = 'sgrna_ids'
    SPECIES = 'species'
    ASSEMBLY = 'assembly'
    MSEQ = 'mseq'
    MSEQ_NO_ADAPT = 'mseq_no_adapt'
    IN_CONST = 'in_const'
    REVC = 'revc'
    OLIGO = 'oligo'
    MUTATION_TYPE = 'mutation_type'
    CDS_PREFIX_LENGTH = 'cds_prefix_length'
    CDS_SUFFIX_LENGTH = 'cds_suffix_length'
    FIRST_CODON_START = 'first_codon_start'
    REF_LEN = 'ref_len'
    ALT_LEN = 'alt_len'
    ALT_REF_DELTA = 'alt_ref_delta'

    def _sql_op(self, op: str) -> str:
        return f"{self.value} {op} ?"

    def ge(self) -> str:
        return self._sql_op('>=')

    def le(self) -> str:
        return self._sql_op('<=')

    def sql_eq(self, value: str | None = None) -> str:
        return f"{self.value}={value or '?'}"


VARIANT_FIELDS = [
    DbFieldName.START,
    DbFieldName.REF,
    DbFieldName.ALT,
    DbFieldName.ID
]


def init_db(conn: Connection) -> None:

    # Load DDL
    with open(get_ddl_path()) as fh:
        ddl = fh.read()

    # Create tables and views
    with cursor(conn) as cur:
        cur.executescript(ddl)


@contextmanager
def get_db_conn() -> Generator[sqlite3.Connection, None, None]:
    with sqlite3.connect(':memory:') as conn:
        yield conn


@contextmanager
def cursor(conn: Connection) -> Generator[Cursor, None, None]:
    cur = conn.cursor()
    try:
        yield cur
    finally:
        cur.close()


def dump_table(conn: Connection, t: DbTableName, fh, sep='\t') -> None:
    for r in conn.execute(f"select * from {t.value}"):
        fh.write(sep.join(str(x) if x is not None else '.' for x in r))
        fh.write('\n')


def dump_all(conn: Connection) -> None:
    for t in DbTableName:
        with open(f"{t.value}_dump.tsv", 'w') as fh:
            dump_table(conn, t, fh)


def get_csv_header(fields: list[DbFieldName]) -> str:
    return ','.join([f.value for f in fields])


def sql_to_str(x: Any) -> str:
    """Serialise Python types expected from a SQL query"""

    if x is None:
        return ''
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    raise TypeError("Unsupported type!")


def select_to_csv(conn: Connection, query: str, field_count: int, header: str, fp: str, sep: str = ',') -> None:
    """
    Stream the results from a SQL query to a CSV file

    Assumptions:
    - the header does not end with a new line character
    """

    n = field_count - 1

    with open(fp, 'w') as fh:

        # Write header
        fh.write(header)
        fh.write('\n')

        with cursor(conn) as cur:
            it = cur.execute(query)

            # Fetch row
            while r := it.fetchone():
                for i in range(n):

                    # Write field
                    fh.write(sql_to_str(r[i]))
                    fh.write(sep)

                # Write last field
                fh.write(sql_to_str(r[n]))
                fh.write('\n')
