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

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd
from .mutated_sequences import BaseMutationCollection
from .oligo_renderer import BaseOligoRenderer
from .pam_protection import PamProtectedReferenceSequence
from .variant import BaseVariant, CustomVariant


@dataclass(init=False)
class CustomVariantOligoRenderer(BaseOligoRenderer):
    __slots__ = {
        'ref_seq',
        'gene_id',
        'transcript_id',
        'adaptor_5',
        'adaptor_3',
        '_oligo_name_prefix'
    }

    def __init__(
        self,
        ref_seq: PamProtectedReferenceSequence,
        gene_id: str,
        transcript_id: str,
        adaptor_5: str,
        adaptor_3: str
    ) -> None:
        super().__init__(ref_seq, gene_id, transcript_id, adaptor_5, adaptor_3)

    def get_oligo_name(self, vcf_alias: str, var_type: int, start: int, ref: Optional[str], alt: Optional[str]) -> str:
        return super()._get_oligo_name(var_type, vcf_alias, start, ref, alt)


@dataclass
class CustomVariantMutation:
    variant: CustomVariant
    sequence: str

    def to_row(self) -> Tuple[Optional[str], Optional[str], int, int, Optional[str], Optional[str], str]:
        var: BaseVariant = self.variant.base_variant
        return (
            self.variant.vcf_alias,
            self.variant.vcf_variant_id,
            var.type.value,
            var.genomic_position.position,
            getattr(var, 'ref', None),
            getattr(var, 'alt', None),
            self.sequence
        )


@dataclass
class CustomVariantMutationCollection(BaseMutationCollection):
    __slots__ = {'df'}

    @classmethod
    def from_variants(cls, mutations: List[CustomVariantMutation]) -> CustomVariantMutationCollection:
        df: pd.DataFrame = pd.DataFrame.from_records([
            mutation.to_row()
            for mutation in mutations
        ], columns=[
            'vcf_alias',
            'vcf_var_id',
            'var_type',
            'mut_position',
            'ref',
            'new',
            'mseq'
        ])

        # Compress table
        df.vcf_alias = df.vcf_alias.astype('category')
        df.vcf_var_id = df.vcf_var_id.astype('string')
        df.var_type = df.var_type.astype(np.int8)
        df.mut_position = df.mut_position.astype(np.int32)
        df.ref = df.ref.astype('category')
        df.new = df.new.astype('category')
        df.mseq = df.mseq.astype('string')

        return cls(df=df)
