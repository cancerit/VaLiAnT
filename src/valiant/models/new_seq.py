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
from typing import List, Optional

from .variant import BaseVariantT, sort_variants


@dataclass(frozen=True)
class VariantGroup:
    __slots__ = ['variants']

    variants: List[BaseVariantT]

    @classmethod
    def from_variants(cls, variants: List[BaseVariantT]):
        return cls(sort_variants(variants))

    def apply(self, start: int, sequence: str, ref_check: bool = False) -> str:
        alt_seq: str = sequence

        for variant in self.variants:
            alt_seq = variant.mutate_from(
                alt_seq, variant.start - start, ref_check=ref_check)

        return alt_seq


@dataclass(frozen=True)
class AltSeqBuilder:
    start: int
    sequence: str
    variant_groups: List[List[BaseVariantT]]
    prefix: str = ''
    suffix: str = ''

    @property
    def ext_sequence(self) -> str:
        return f"{self.prefix}{self.sequence}{self.suffix}"

    def get_alt(self, extend: bool = False, variant_layer: Optional[int] = None) -> str:
        n: int = len(self.variant_groups)
        last_group_index: int = (
            variant_layer if variant_layer is not None else
            n - 1
        )
        if last_group_index >= n:
            raise ValueError(f"Invalid variant group index {variant_layer}!")

        alt_seq: str = self.ext_sequence if extend else self.sequence
        for g in self.variant_groups[:last_group_index + 1]:
            alt_seq = g.apply(self.start, alt_seq)

        return alt_seq

    def get_alt_length(self, extend: bool = False, variant_layer: Optional[int] = None) -> int:
        return len(self.get_alt(extend=extend, variant_layer=variant_layer))
