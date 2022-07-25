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
from typing import List, Optional, Tuple, FrozenSet
import numpy as np
import pandas as pd
from ..constants import META_PAM_MUT_SGRNA_ID, META_VCF_VAR_IN_CONST
from ..sgrna_utils import sgrna_ids_to_string
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


@dataclass(frozen=True)
class CustomVariantMutation:
    __slots__ = ['variant', 'sequence', 'overlaps_constant_region', 'sgrna_ids']

    variant: CustomVariant
    sequence: str
    overlaps_constant_region: bool
    sgrna_ids: FrozenSet[str]

    @property
    def ref(self) -> Optional[str]:
        return getattr(self.variant.base_variant, 'ref', None)

    @property
    def alt(self) -> Optional[str]:
        return getattr(self.variant.base_variant, 'alt', None)

    @property
    def ref_length(self) -> int:
        return len(self.ref) if self.ref else 0

    def to_row(self) -> Tuple[Optional[str], Optional[str], int, int, Optional[str], Optional[str], str, str]:
        var: BaseVariant = self.variant.base_variant
        return (
            self.variant.vcf_alias,
            self.variant.vcf_variant_id,
            var.type.value,
            var.genomic_position.position,
            self.ref,
            self.alt,
            self.sequence,
            sgrna_ids_to_string(self.sgrna_ids) if self.sgrna_ids else '',
            1 if self.overlaps_constant_region else 0
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
            'mseq',
            META_PAM_MUT_SGRNA_ID,
            META_VCF_VAR_IN_CONST
        ])

        # Compress table
        df.vcf_alias = df.vcf_alias.astype('category')
        df.vcf_var_id = df.vcf_var_id.astype('string')
        df.var_type = df.var_type.astype(np.int8)
        df.mut_position = df.mut_position.astype(np.int32)
        df.ref = df.ref.astype('category')
        df.new = df.new.astype('category')
        df.mseq = df.mseq.astype('string')
        df[META_PAM_MUT_SGRNA_ID] = df[META_PAM_MUT_SGRNA_ID].astype('string')
        df[META_VCF_VAR_IN_CONST] = df[META_VCF_VAR_IN_CONST].astype(pd.Int8Dtype())

        return cls(df=df)
