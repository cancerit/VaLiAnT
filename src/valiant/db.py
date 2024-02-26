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
    POS_R = 'pos_r'
    POS_A = 'pos_a'
    REF_A = 'ref_a'
    ALT_A = 'alt_a'
    MUTATOR = 'mutator'
    CODON_REF_A = 'codon_ref_a'
    CODON_ALT_A = 'codon_alt_a'
    AA_REF = 'aa_ref'
    AA_ALT = 'aa_alt'
    REF_START = 'ref_start'
    SGRNA_IDS = 'sgrna_ids'
    IN_CONST = 'in_const'
    OLIGO = 'oligo'
    MUTATION_TYPE = 'mutation_type'
    CDS_PREFIX_LENGTH = 'cds_prefix_length'
    CDS_SUFFIX_LENGTH = 'cds_suffix_length'
    FIRST_CODON_START = 'first_codon_start'
    REF_LEN = 'ref_len'
    ALT_LEN = 'alt_len'

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
