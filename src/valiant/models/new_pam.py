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

from __future__ import annotations

import abc
from dataclasses import dataclass, field
import logging
from typing import Dict, FrozenSet, Generic, List, Optional, TypeVar

from .alt_seq_builder import AltSeqBuilder, AltSeqBuilderT
from .base import GenomicRange, StrandedPositionRange
from .cds_alt_seq_builder import CdsAltSeqBuilder
from .codon_table import CodonTable
from .dna_str import DnaStr
from .pam_protection import PamVariant
from .sequences import ReferenceSequence
from .variant import BaseVariant, BaseVariantT
from .variant_group import VariantGroup
from ..enums import MutationType

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
    pam_seq: str = field(init=False)
    ext_pam_seq: str = field(init=False)

    def __post_init__(self) -> None:

        # Cache the ALT sequence
        for attr, ext in [
            ('pam_seq', False),
            ('ext_pam_seq', True)
        ]:
            object.__setattr__(self, attr, self.get_pam_seq(
                extend=ext, ref_check=True))

    @property
    def pos_range(self) -> GenomicRange:
        return self.ab.gr

    @property
    def start(self) -> int:
        return self.ab.start

    @property
    def ref_seq(self) -> str:
        return self.ab.ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.ext_pam_seq

    @property
    def alt_seq(self) -> str:
        return self.pam_seq

    @property
    def bg_variants(self) -> List[BaseVariant]:
        return self.ab.variant_groups[LAYER_BG].variants

    @property
    def pam_variants(self) -> List[PamVariant]:
        return self.ab.variant_groups[LAYER_PAM].variants

    def get_pam_variants_in_range(self, spr: StrandedPositionRange) -> List[PamVariant]:
        return self.ab.get_variants(LAYER_PAM, genomic_range=spr)

    def get_sgrna_ids_in_range(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return frozenset({
            variant.sgrna_id
            for variant in self.get_pam_variants_in_range(spr)
        })

    def _get_alt(self, layer: int, extend: bool = False, ref_check: bool = False) -> str:
        return self.ab.get_alt(variant_layer=layer, extend=extend, ref_check=ref_check)

    def get_bg_seq(self, extend: bool = False, ref_check: bool = False) -> str:
        return self._get_alt(LAYER_BG, extend=extend, ref_check=ref_check)

    def get_pam_seq(self, extend: bool = False, ref_check: bool = False) -> str:
        return self._get_alt(LAYER_PAM, extend=extend, ref_check=ref_check)

    def mutate(self, variant: BaseVariant) -> str:
        return self.ab.mutate_alt(variant)

    def get_variant_corrected_ref(self, variant: BaseVariant) -> Optional[str]:
        return variant.get_corrected_ref(self.pam_seq, self.start)


@dataclass(frozen=True)
class PamBgAltSeqBuilder(BasePamBgAltSeqBuilder[AltSeqBuilder]):

    @classmethod
    def from_ref(
        cls,
        ref_range: GenomicRange,
        ref_seq: str,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant]
    ) -> PamBgAltSeqBuilder:
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
    ) -> PamBgAltSeqBuilder:
        return cls(AltSeqBuilder(
            ref_seq.genomic_range,
            DnaStr(ref_seq.sequence),
            _get_variant_groups(bg_variants, pam_variants)))

    def get_pam_sub(self, gr: GenomicRange) -> PamBgAltSeqBuilder:
        return PamBgAltSeqBuilder(self.ab.get_sub(gr))


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
    ) -> CdsPamBgAltSeqBuilder:
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
    ) -> CdsPamBgAltSeqBuilder:
        return cls.from_ref(
            ref_seq.genomic_range,
            ref_seq.sequence,
            bg_variants,
            pam_variants,
            cds_prefix=cds_prefix,
            cds_suffix=cds_suffix)

    @classmethod
    def from_noncds(cls, x: PamBgAltSeqBuilder, prefix: str, suffix: str) -> CdsPamBgAltSeqBuilder:
        return cls(CdsAltSeqBuilder.from_builder(x.ab, prefix, suffix))

    @property
    def frame(self) -> int:
        return self.ab.frame

    @property
    def cds_prefix(self) -> str:
        return self.ab.cds_prefix

    @property
    def cds_suffix(self) -> str:
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

    def overlaps_bg(self, variant: BaseVariant) -> bool:
        """Does the variant overlap with any of the background variants?"""

        return self.ab.overlaps_layer(LAYER_BG, variant)

    def get_codon_indices_in_range(self, spr: StrandedPositionRange) -> List[int]:
        return self.ab.get_codon_indices_in_range(spr)

    def any_non_syn_background(self, codon_table: CodonTable) -> bool:
        return any(
            self.ab.is_variant_nonsynonymous(codon_table, variant)
            for variant in self.bg_variants
        )

    def any_frame_shifting_background(self) -> bool:
        return any(
            variant.is_frame_shifting
            for variant in self.bg_variants
        )

    def is_background_valid(self, codon_table: CodonTable, allow_frame_shift: bool, allow_non_syn: bool) -> bool:
        any_fs: bool = False
        any_ns: bool = False
        for variant in self.bg_variants:
            if self.ab.is_variant_nonsynonymous(codon_table, variant):
                fs = variant.is_frame_shifting
                any_ns = True
                if fs:
                    any_fs = True
                    prefix = "Frame shifting"
                else:
                    prefix = "Non-synonymous"
                logging.warning(f"{prefix} background variant at {variant.genomic_position}!")

        return not (
            (not allow_non_syn and any_ns) or
            (not allow_frame_shift and any_fs)
        )
