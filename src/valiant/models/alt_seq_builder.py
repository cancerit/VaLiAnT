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

from dataclasses import dataclass
from typing import List, Optional, TypeVar

from .base import GenomicPosition, GenomicRange
from .dna_str import DnaStr
from .uint_range import UIntRange
from .variant import BaseVariant, BaseVariantT as VariantT
from .variant_group import VariantGroup
from ..utils import has_duplicates

RangeT = TypeVar('RangeT', bound='StrandedPositionRange')


class InvalidVariantGroupIndex(IndexError):
    pass


def is_valid_index(range_length: int, value: int) -> bool:
    return 0 <= value < range_length


@dataclass(frozen=True)
class AltSeqBuilder:
    __slots__ = ['gr', 'ref_seq', 'variant_groups']

    gr: GenomicRange
    ref_seq: DnaStr
    variant_groups: List[VariantGroup]

    # Genomic position properties

    @property
    def start(self) -> int:
        return self.gr.start

    @property
    def end(self) -> int:
        return self.gr.end

    # Reference sequence properties

    @property
    def seq_length(self) -> int:
        return len(self.ref_seq)

    @property
    def ext_seq_length(self) -> int:
        return self.seq_length

    @property
    def ext_ref_seq(self) -> str:
        return self.ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.get_alt(extend=True)

    # Variant group properties

    @property
    def variant_group_count(self) -> int:
        return len(self.variant_groups)

    def _extend(self, s: str) -> str:
        return s

    def is_variant_group_index_valid(self, variant_group_index: int) -> bool:
        return is_valid_index(self.variant_group_count, variant_group_index)

    def validate_variant_group_index(self, variant_group_index: int) -> None:
        if not self.is_variant_group_index_valid(variant_group_index):
            raise InvalidVariantGroupIndex(f"Invalid variant group index {variant_group_index}!")

    def are_variant_group_indices_valid(self, variant_group_indices: List[int]) -> bool:
        return all(map(self.is_variant_group_index_valid, variant_group_indices))

    def validate_variant_group_indices(self, variant_group_indices: List[int]) -> None:
        if has_duplicates(variant_group_indices):
            raise ValueError("Duplicate variant group indices!")
        if not self.are_variant_group_indices_valid(variant_group_indices):
            raise InvalidVariantGroupIndex("Invalid variant group indices!")

    def get_variants(self, variant_group_index: int, genomic_range: Optional[RangeT] = None) -> List[VariantT]:
        def ft(variant: VariantT) -> bool:
            return variant.in_range(genomic_range)  # type: ignore

        variants = self.get_variant_group(variant_group_index).variants
        return list(filter(ft, variants)) if genomic_range else variants

    def get_variant_positions(self, variant_group_index: int) -> List[GenomicPosition]:
        return [
            variant.genomic_position
            for variant in self.get_variants(variant_group_index)
        ]

    def get_variant_ref_ranges(self, variant_group_index: int) -> List[GenomicRange]:
        return [
            variant.ref_range
            for variant in self.get_variants(variant_group_index)
        ]

    def __post_init__(self) -> None:
        if self.seq_length != len(self.gr):
            raise ValueError("Mismatching position range and sequence length!")

    def get_alt(
        self,
        extend: bool = False,
        variant_layer: Optional[int] = None,
        ref_check: bool = False,
        correct_alt: bool = False
    ) -> str:
        last_group_index: int = (
            variant_layer if variant_layer is not None else
            self.variant_group_count - 1
        )
        self.validate_variant_group_index(last_group_index)

        alt_seq: str = self.ref_seq
        alt_offset: int = 0

        # Only allow reference checks on the first layer
        # TODO: should warnings be raised for failed reference checks on the following layers?
        for i, g in enumerate(self.variant_groups[:last_group_index + 1]):
            alt_seq, layer_alt_offset = g.apply(
                self.start,
                alt_seq,
                alt_offset=alt_offset,
                ref_check=ref_check if i == 0 else False)
            if correct_alt:
                alt_offset += layer_alt_offset

        return self._extend(alt_seq) if extend else alt_seq

    def get_alt_length(self, extend: bool = False, variant_layer: Optional[int] = None) -> int:
        return len(self.get_alt(extend=extend, variant_layer=variant_layer))

    def _get_sub_seq(self, r: UIntRange) -> DnaStr:
        nr: UIntRange = r - self.start
        return DnaStr(self.ref_seq[nr.to_slice()])

    def _get_sub_groups(self, r: UIntRange) -> List[VariantGroup]:
        return [
            g.get_sub(r)
            for g in self.variant_groups
        ]

    def get_sub(self, gr: GenomicRange) -> AltSeqBuilder:
        r = gr.to_uintr()
        return AltSeqBuilder(
            gr,
            self._get_sub_seq(r),
            self._get_sub_groups(r))

    def get_variant_group(self, variant_group_index: int) -> VariantGroup:
        self.validate_variant_group_index(variant_group_index)
        return self.variant_groups[variant_group_index]

    def overlaps_layer(self, variant_group_index: int, variant: BaseVariant) -> bool:
        return self.get_variant_group(variant_group_index).overlaps(variant)

    def mutate_alt(self, variant: BaseVariant, variant_layer: Optional[int] = None, ref_check: bool = False) -> str:
        if not variant.in_range(self.gr):
            raise RuntimeError("Variant out of bounds!")
        # TODO: verify support for variants that begin before the reference start
        return variant.mutate(
            self.get_alt(variant_layer=variant_layer, ref_check=False),
            self.start,
            ref_check=ref_check)


AltSeqBuilderT = TypeVar('AltSeqBuilderT', bound=AltSeqBuilder)
