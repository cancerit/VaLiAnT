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
from typing import List, Optional, TypeVar

from ..utils import has_duplicates
from .base import GenomicPosition, GenomicRange, StrandedPositionRange
from .dna_str import DnaStr
from .uint_range import UIntRange
from .variant import BaseVariantT as VariantT
from .variant_group import VariantGroup


RangeT = TypeVar('RangeT', bound='StrandedPositionRange')


class InvalidVariantGroupIndex(IndexError):
    pass


def is_valid_index(range_length: int, value: int) -> bool:
    return 0 <= value < range_length


@dataclass(frozen=True)
class AltSeqBuilder:
    gr: GenomicRange
    sequence: DnaStr
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
    def ext_sequence(self) -> str:
        return self.sequence

    @property
    def seq_length(self) -> int:
        return len(self.sequence)

    @property
    def ext_seq_length(self) -> int:
        return self.seq_length

    @property
    def ext_ref_seq(self) -> str:
        return self.ext_sequence

    @property
    def ext_alt_seq(self) -> str:
        return self.get_alt(extend=True)

    # Variant group properties

    @property
    def variant_group_count(self) -> int:
        return len(self.variant_groups)

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
            return (
                variant.ref_start_pos.in_range(genomic_range) or
                variant.ref_end_pos.in_range(genomic_range)
            )

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

    def get_alt(self, extend: bool = False, variant_layer: Optional[int] = None) -> str:
        last_group_index: int = (
            variant_layer if variant_layer is not None else
            self.variant_group_count - 1
        )
        self.validate_variant_group_index(last_group_index)

        alt_seq: str = self.ext_sequence if extend else self.sequence
        for g in self.variant_groups[:last_group_index + 1]:
            alt_seq = g.apply(self.start, alt_seq)

        return alt_seq

    def get_alt_length(self, extend: bool = False, variant_layer: Optional[int] = None) -> int:
        return len(self.get_alt(extend=extend, variant_layer=variant_layer))

    def get_sub(self, r: UIntRange) -> 'AltSeqBuilder':
        nr: UIntRange = r - self.start

        return AltSeqBuilder(
            r.start,
            self.sequence[nr.to_slice()],
            [
                g.get_sub(r)
                for g in self.variant_groups
            ])

    def get_variant_group(self, variant_group_index: int) -> VariantGroup:
        self.validate_variant_group_index(variant_group_index)
        return self.variant_groups[variant_group_index]
