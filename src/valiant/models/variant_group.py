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

from .uint_range import UIntRange
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

    def get_sub(self, r: UIntRange) -> 'VariantGroup':
        def var_in_range(variant: BaseVariantT) -> bool:
            if variant.start in r:
                return True
            if variant.ref_length > 1:
                var_end = variant.start + variant.ref_length - 1
                if var_end in r:
                    return True
            return False

        return VariantGroup(list(filter(var_in_range, self.variants)))
