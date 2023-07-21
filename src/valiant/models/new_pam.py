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

from .alt_seq_builder import AltSeqBuilder
from .cds_alt_seq_builder import CdsAltSeqBuilder
from .dna_str import DnaStr
from .pam_protection import PamVariant
from .sequences import ReferenceSequence
from .variant import BaseVariantT
from .variant_group import VariantGroup


LAYER_BG = 0
LAYER_PAM = 1


@dataclass(frozen=True, init=False)
class PamBgAltSeqBuilder(AltSeqBuilder):
    @classmethod
    def from_ref_seq(
        cls,
        ref_seq: ReferenceSequence,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant]
    ) -> 'PamBgAltSeqBuilder':
        # Beware: maintain the variant layer order!
        return cls(
            ref_seq.genomic_range,
            DnaStr(ref_seq.sequence), [
                VariantGroup.from_variants(variants)
                for variants in [
                    bg_variants,  # LAYER_BG
                    pam_variants  # LAYER_PAM
                ]
            ])

    @property
    def bg_variants(self) -> List[BaseVariantT]:
        return self.variant_groups[LAYER_BG].variants

    @property
    def pam_variants(self) -> List[PamVariant]:
        return self.variant_groups[LAYER_PAM].variants

    def get_bg_seq(self, extend: bool = False) -> str:
        return self.get_alt(extend=extend, variant_layer=LAYER_BG)

    def get_pam_seq(self, extend: bool = False) -> str:
        return self.get_alt(extend=extend, variant_layer=LAYER_PAM)


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
            ref_seq.genomic_range,
            DnaStr(ref_seq.sequence), [
                VariantGroup.from_variants(variants)
                for variants in [
                    bg_variants,
                    pam_variants
                ]
            ],
            prefix=DnaStr(cds_prefix),
            suffix=DnaStr(cds_suffix))

    @property
    def pam_bg_codon_clash(self) -> bool:
        return self.variant_group_codon_clash([
            LAYER_BG,
            LAYER_PAM
        ])
