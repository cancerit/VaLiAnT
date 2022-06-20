########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from typing import Callable, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from .options import Options
from .pam_protection import PamProtectedReferenceSequence
from ..constants import REVCOMP_OLIGO_NAME_SUFFIX
from ..enums import VariantType
from ..utils import get_constant_category, reverse_complement, get_source_type_column


var_type_sub: int = VariantType.SUBSTITUTION.value
var_type_del: int = VariantType.DELETION.value
var_type_ins: int = VariantType.INSERTION.value


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
            raise ValueError("Invalid insertion: missing alternative!")
        return f"{oligo_name_prefix}{start}_{alt}_{source}"

    else:
        if not ref:
            if var_type == var_type_sub:
                raise ValueError(f"Invalid substitution: missing reference!")
            if var_type == var_type_del:
                raise ValueError(f"Invalid deletion: missing reference!")
            else:
                raise ValueError("Invalid variant type!")

        ref_len: int = len(ref)
        end: int = (start + ref_len - 1) if ref_len > 1 else start

        # Substitution
        if var_type == var_type_sub:
            if not alt:
                raise ValueError("Invalid substitution: missing alternative!")
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

        raise ValueError("Invalid variant type!")


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

    ref_seq: PamProtectedReferenceSequence
    gene_id: str
    transcript_id: str
    adaptor_5: str
    adaptor_3: str
    _oligo_name_prefix: str

    def __init__(
        self,
        ref_seq: PamProtectedReferenceSequence,
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
        return self.ref_seq.genomic_range.chromosome

    @property
    def start(self) -> int:
        return self.ref_seq.genomic_range.start

    @property
    def strand(self) -> str:
        return self.ref_seq.genomic_range.strand

    @property
    def _constant_fields(self) -> List[Tuple[str, Any]]:
        return [
            ('gene_id', self.gene_id),
            ('transcript_id', self.transcript_id),
            ('ref_chr', self.chromosome),
            ('ref_strand', self.ref_seq.genomic_range.strand),
            ('ref_seq', self.ref_seq.sequence),
            ('pam_seq', self.ref_seq.pam_protected_sequence)
        ]

    def _render_mutated_sequence(self, mseq: str) -> str:
        return f"{self.adaptor_5}{mseq}{self.adaptor_3}"

    def _render_mutated_sequence_rc(self, mseq: str) -> str:
        return f"{self.adaptor_5}{reverse_complement(mseq)}{self.adaptor_3}"

    def _get_renderer(self, reverse_complement: bool) -> Callable[[str], str]:
        if reverse_complement:
            if self.strand != '-':
                raise ValueError("Attempted reverse complement of plus strand sequence!")
            return self._render_mutated_sequence_rc
        else:
            return self._render_mutated_sequence

    def _get_oligo_name(self, var_type: int, source: str, start: int, ref: Optional[str], alt: Optional[str]) -> str:
        return get_oligo_name(self._oligo_name_prefix, var_type, source, start, ref, alt)

    def get_metadata_table(self, df: pd.DataFrame, options: Options) -> pd.DataFrame:
        if set(df.columns.array) < {'oligo_name', 'mut_position', 'ref', 'new', 'mutator', 'mseq'}:
            raise ValueError("Invalid mutation metadata data frame!")

        rown: int = df.shape[0]

        # Compress mutation metadata
        df.oligo_name = df.oligo_name.astype('string')
        df.mutator = df.mutator.astype('category')

        # Add global metadata
        df['ref_start'] = np.int32(self.ref_seq.genomic_range.start)
        df['ref_end'] = np.int32(self.ref_seq.genomic_range.end)
        for col_name, col_value in self._constant_fields:
            df[col_name] = get_constant_category(col_value, rown)

        # Add reverse complement information to the metadata
        rc: bool = self.ref_seq.genomic_range.strand == '-' and options.revcomp_minus_strand
        df['revc'] = np.repeat(rc, rown).astype(np.int8)
        if rc:
            df.oligo_name = df.oligo_name + REVCOMP_OLIGO_NAME_SUFFIX

        # Render full oligonucleotide sequences
        df.mseq = df.mseq.apply(self._get_renderer(rc)).astype('string')

        # Set sequence source type
        df['src_type'] = get_source_type_column('ref', rown)

        return df
