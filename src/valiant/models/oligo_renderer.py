########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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

from dataclasses import dataclass
from functools import partial
from typing import List, NoReturn, Optional, Set, Tuple, Any
import numpy as np
import pandas as pd
from .new_pam import PamBgAltSeqBuilderT
from .options import Options
from ..constants import META_MSEQ, META_MSEQ_NO_ADAPT, META_MSEQ_NO_ADAPT_NO_RC, META_MUT_POSITION, META_MUTATOR, META_NEW, META_OLIGO_NAME, META_REF, META_REF_END, META_REF_START, REVCOMP_OLIGO_NAME_SUFFIX
from ..enums import VariantType
from ..utils import get_constant_category, reverse_complement, get_source_type_column, validate_strand


ERR_MISSING_ALT = 'missing alternative'
ERR_MISSING_REF = 'missing reference'

var_type_sub: int = VariantType.SUBSTITUTION.value
var_type_del: int = VariantType.DELETION.value
var_type_ins: int = VariantType.INSERTION.value

MIN_METADATA: Set[str] = {
    META_OLIGO_NAME,
    META_MUT_POSITION,
    META_REF,
    META_NEW,
    META_MUTATOR,
    META_MSEQ
}


def raise_err_var_type(var_type_value: int, msg: str) -> NoReturn:
    raise ValueError(f"Invalid {VariantType(var_type_value).name.lower()}: {msg}!")


def raise_err_invalid_var_type() -> NoReturn:
    raise ValueError("Invalid variant type!")


def get_oligo_name(
    oligo_name_prefix: str,
    var_type: int,
    source: str,
    start: int,
    ref: Optional[str],
    alt: Optional[str]
) -> str:

    # Insertion
    if var_type == var_type_ins:
        if not alt:
            raise_err_var_type(var_type, ERR_MISSING_ALT)
        return f"{oligo_name_prefix}{start}_{alt}_{source}"

    else:
        if not ref:
            if var_type == var_type_sub:
                raise_err_var_type(var_type, ERR_MISSING_REF)
            if var_type == var_type_del:
                raise_err_var_type(var_type, ERR_MISSING_REF)
            else:
                raise_err_invalid_var_type()

        ref_len: int = len(ref)
        end: int = (start + ref_len - 1) if ref_len > 1 else start

        # Substitution
        if var_type == var_type_sub:
            if not alt:
                raise_err_var_type(var_type, ERR_MISSING_ALT)
            return (
                f"{oligo_name_prefix}{start}_{ref}>{alt}_{source}" if end == start else
                f"{oligo_name_prefix}{start}_{end}_{ref}>{alt}_{source}"
            )

        # Deletion
        if var_type == var_type_del:
            return (
                f"{oligo_name_prefix}{start}_{source}" if end == start else
                f"{oligo_name_prefix}{start}_{end}_{source}"
            )

        raise_err_invalid_var_type()


def get_oligo_name_from_row(oligo_name_prefix: str, r) -> str:
    return get_oligo_name(
        oligo_name_prefix,
        r.var_type,
        r.mutator,
        r.mut_position,
        r.ref if not pd.isnull(r.ref) else None,
        r.new if not pd.isnull(r.new) else None)


def get_oligo_names_from_dataframe(oligo_name_prefix: str, df: pd.DataFrame) -> pd.Series:
    f = partial(get_oligo_name_from_row, oligo_name_prefix)
    return pd.Series(df.apply(f, axis=1), dtype='string')


@dataclass(init=False)
class BaseOligoRenderer:
    __slots__ = {
        'ref_seq',
        'gene_id',
        'transcript_id',
        'adaptor_5',
        'adaptor_3',
        '_oligo_name_prefix'
    }

    ref_seq: PamBgAltSeqBuilderT
    gene_id: str
    transcript_id: str
    adaptor_5: str
    adaptor_3: str
    _oligo_name_prefix: str

    def __init__(
        self,
        ref_seq: PamBgAltSeqBuilderT,
        gene_id: str,
        transcript_id: str,
        adaptor_5: str,
        adaptor_3: str
    ) -> None:
        self.ref_seq = ref_seq
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.adaptor_5 = adaptor_5
        self.adaptor_3 = adaptor_3

        self._oligo_name_prefix = (
            f"{self.transcript_id}.{self.gene_id}" if self.transcript_id and self.gene_id else
            'NO_TRANSCRIPT'
        ) + f"_{self.chromosome}:"

    @property
    def chromosome(self) -> str:
        return self.ref_seq.pos_range.chromosome

    @property
    def start(self) -> int:
        return self.ref_seq.pos_range.start

    @property
    def strand(self) -> str:
        return self.ref_seq.pos_range.strand

    @property
    def has_adaptors(self) -> bool:
        return bool(self.adaptor_5) or bool(self.adaptor_3)

    @property
    def _constant_fields(self) -> List[Tuple[str, Any]]:
        # TODO: background variants in `ref_seq`?
        return [
            ('gene_id', self.gene_id),
            ('transcript_id', self.transcript_id),
            ('ref_chr', self.chromosome),
            ('ref_strand', self.ref_seq.pos_range.strand),
            ('ref_seq', self.ref_seq.ref_seq),
            ('pam_seq', self.ref_seq.get_pam_seq())
        ]

    def _add_adaptors(self, mseq: str) -> str:
        return f"{self.adaptor_5}{mseq}{self.adaptor_3}"

    def _should_apply_reverse_complement(self, apply_reverse_complement: bool) -> bool:
        validate_strand(self.strand)

        if apply_reverse_complement and self.strand != '-':
            raise ValueError("Attempted reverse complement of plus strand sequence!")

        return apply_reverse_complement

    def _get_oligo_name(self, var_type: int, source: str, start: int, ref: Optional[str], alt: Optional[str]) -> str:
        return get_oligo_name(self._oligo_name_prefix, var_type, source, start, ref, alt)

    def get_oligo_names_from_dataframe(self, df: pd.DataFrame) -> pd.Series:
        return get_oligo_names_from_dataframe(self._oligo_name_prefix, df)

    def _get_mutated_sequence(self, mseq: str) -> str:
        return mseq

    def _render_mutated_sequence_no_adaptors(self, mseq: str) -> str:
        return self._get_mutated_sequence(mseq)

    def get_metadata_table(self, df: pd.DataFrame, options: Options) -> pd.DataFrame:
        if set(df.columns.array) < MIN_METADATA:
            raise ValueError("Invalid mutation metadata data frame!")

        rown: int = df.shape[0]

        # Compress mutation metadata
        df.oligo_name = df.oligo_name.astype('string')
        df.mutator = df.mutator.astype('category')

        # Add global metadata
        df[META_REF_START] = np.int32(self.ref_seq.pos_range.start)
        df[META_REF_END] = np.int32(self.ref_seq.pos_range.end)
        for col_name, col_value in self._constant_fields:
            df[col_name] = get_constant_category(col_value, rown)

        # Add reverse complement information to the metadata
        rc: bool = self._should_apply_reverse_complement(options.revcomp_minus_strand)
        df['revc'] = np.repeat(rc, rown).astype(np.int8)
        if rc:
            df.oligo_name = df.oligo_name + REVCOMP_OLIGO_NAME_SUFFIX

        # Render oligonucleotide sequences without adaptors (never reverse complement)
        df[META_MSEQ_NO_ADAPT_NO_RC] = df[META_MSEQ].apply(
            self._render_mutated_sequence_no_adaptors).astype('string')

        # Render oligonucleotide sequences without adaptors
        df[META_MSEQ_NO_ADAPT] = (
            df[META_MSEQ_NO_ADAPT_NO_RC].apply(reverse_complement).astype('string') if rc else
            df[META_MSEQ_NO_ADAPT_NO_RC]
        )

        # Render full oligonucleotide sequences
        df[META_MSEQ] = (
            df[META_MSEQ_NO_ADAPT].apply(self._add_adaptors).astype('string') if self.has_adaptors else
            df[META_MSEQ_NO_ADAPT]
        )

        # Set sequence source type
        df['src_type'] = get_source_type_column('ref', rown)

        return df
