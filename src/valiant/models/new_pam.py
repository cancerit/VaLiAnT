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

from dataclasses import dataclass
from typing import List

from .new_seq import AltSeqBuilder, CdsAltSeqBuilder, DnaStr, VariantGroup
from .pam_protection import PamVariant
from .sequences import ReferenceSequence
from .variant import BaseVariantT


@dataclass(frozen=True, init=False)
class PamBgAltSeqBuilder(AltSeqBuilder):
    @classmethod
    def from_ref_seq(
        cls,
        ref_seq: ReferenceSequence,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant]
    ) -> 'PamBgAltSeqBuilder':
        return cls(
            ref_seq.start,
            DnaStr(ref_seq.sequence), [
                VariantGroup.from_variants(variants)
                for variants in [
                    bg_variants,
                    pam_variants
                ]
            ])


@dataclass(frozen=True, init=False)
class CdsPamBgAltSeqBuilder(CdsAltSeqBuilder):
    @classmethod
    def from_ref_seq(
        cls,
        ref_seq: ReferenceSequence,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant],
        cds_prefix: str = '',
        cds_suffix: str = ''
    ) -> 'CdsPamBgAltSeqBuilder':
        return cls(
            ref_seq.start,
            DnaStr(ref_seq.sequence), [
                VariantGroup.from_variants(variants)
                for variants in [
                    bg_variants,
                    pam_variants
                ]
            ],
            prefix=DnaStr(cds_prefix),
            suffix=DnaStr(cds_suffix))
