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

import abc
from dataclasses import dataclass
from typing import Dict, Generic, List, TypeVar

from ..enums import MutationType

from .alt_seq_builder import AltSeqBuilderT, AltSeqBuilder
from .base import GenomicRange, StrandedPositionRange
from .cds_alt_seq_builder import CdsAltSeqBuilder
from .codon_table import CodonTable
from .dna_str import DnaStr
from .pam_protection import PamVariant
from .sequences import ReferenceSequence
from .variant import BaseVariantT
from .variant_group import VariantGroup


LAYER_BG = 0
LAYER_PAM = 1
LAYERS = [
    LAYER_BG,
    LAYER_PAM
]


def _get_variant_groups(
    bg_variants: List[BaseVariantT],
    pam_variants: List[PamVariant]
) -> List[VariantGroup]:
    # Beware: maintain the variant layer order!
    return [
        VariantGroup.from_variants(variants)
        for variants in [
            bg_variants,  # LAYER_BG
            pam_variants  # LAYER_PAM
        ]
    ]


PamBgAltSeqBuilderT = TypeVar('PamBgAltSeqBuilderT', bound='BasePamBgAltSeqBuilder')


@dataclass(frozen=True)
class BasePamBgAltSeqBuilder(abc.ABC, Generic[AltSeqBuilderT]):
    ab: AltSeqBuilderT

    @property
    def pos_range(self) -> StrandedPositionRange:
        return self.ab.gr

    @property
    def start(self) -> int:
        return self.ab.start

    @property
    def ref_seq(self) -> str:
        return self.ab.ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.ab.ext_alt_seq

    @property
    def alt_seq(self) -> str:
        return self.ab.get_alt(extend=False)

    @property
    def bg_variants(self) -> List[BaseVariantT]:
        return self.ab.variant_groups[LAYER_BG].variants

    @property
    def pam_variants(self) -> List[PamVariant]:
        return self.ab.variant_groups[LAYER_PAM].variants

    def get_pam_variants_in_range(self, spr: StrandedPositionRange) -> List[PamVariant]:
        return self.ab.get_variants(LAYER_PAM, genomic_range=spr)

    def get_bg_seq(self, extend: bool = False, ref_check: bool = False) -> str:
        return self.ab.get_alt(extend=extend, variant_layer=LAYER_BG, ref_check=ref_check)

    def get_pam_seq(self, extend: bool = False, ref_check: bool = False) -> str:
        return self.ab.get_alt(extend=extend, variant_layer=LAYER_PAM, ref_check=ref_check)


@dataclass(frozen=True)
class PamBgAltSeqBuilder(BasePamBgAltSeqBuilder[AltSeqBuilder]):

    @classmethod
    def from_ref(
        cls,
        ref_range: GenomicRange,
        ref_seq: str,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant]
    ) -> 'PamBgAltSeqBuilder':
        return cls(AltSeqBuilder(
            ref_range,
            DnaStr(ref_seq),
            _get_variant_groups(bg_variants, pam_variants)))

    @classmethod
    def from_ref_seq(
        cls,
        ref_seq: ReferenceSequence,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant]
    ) -> 'PamBgAltSeqBuilder':
        return cls(AltSeqBuilder(
            ref_seq.genomic_range,
            ref_seq.sequence,
            _get_variant_groups(bg_variants, pam_variants)))


@dataclass(frozen=True)
class CdsPamBgAltSeqBuilder(BasePamBgAltSeqBuilder[CdsAltSeqBuilder]):

    @classmethod
    def from_ref(
        cls,
        ref_range: GenomicRange,
        ref_seq: str,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant],
        cds_prefix: str = '',
        cds_suffix: str = ''
    ) -> 'CdsPamBgAltSeqBuilder':
        return cls(CdsAltSeqBuilder(
            ref_range,
            DnaStr(ref_seq),
            _get_variant_groups(bg_variants, pam_variants),
            DnaStr(cds_prefix),
            DnaStr(cds_suffix)))

    @classmethod
    def from_ref_seq(
        cls,
        ref_seq: ReferenceSequence,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant],
        cds_prefix: str = '',
        cds_suffix: str = ''
    ) -> 'CdsPamBgAltSeqBuilder':
        return cls.from_ref(
            ref_seq.genomic_range,
            ref_seq.sequence,
            bg_variants,
            pam_variants,
            cds_prefix=cds_prefix,
            cds_suffix=cds_suffix)

    @property
    def frame(self) -> int:
        return self.ab.frame

    @property
    def cds_prefix(self) -> int:
        return self.ab.cds_prefix

    @property
    def cds_suffix(self) -> int:
        return self.ab.cds_suffix

    @property
    def cds_prefix_length(self) -> int:
        return self.ab.cds_prefix_length

    @property
    def cds_suffix_length(self) -> int:
        return self.ab.cds_suffix_length

    @property
    def pam_bg_codon_clash(self) -> bool:
        return self.ab.variant_group_codon_clash(LAYERS)

    @property
    def contains_same_codon_variants(self) -> bool:
        return self.ab.contains_same_codon_variants(LAYER_PAM)

    @property
    def contains_same_codon_variants(self) -> bool:
        return self.ab.contains_same_codon_variants(LAYER_PAM)

    def get_pam_variant_codon_indices(self) -> List[int]:
        return self.ab.get_variant_group_codon_indices(LAYER_PAM)

    def get_indexed_alt_codons(self) -> Dict[int, str]:
        return self.ab.get_indexed_alt_codons()

    def log_same_codon_variants(self) -> None:
        return self.ab.log_same_codon_variants(LAYER_PAM)

    def get_variant_mutation_types(self, codon_table: CodonTable, **kwargs) -> List[MutationType]:
        return self.ab.get_variant_mutation_types(codon_table, LAYER_PAM, **kwargs)

    def get_ref_codon_index(self, position: int) -> int:
        return self.ab.get_ref_codon_index(position)
