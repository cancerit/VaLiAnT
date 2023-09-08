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

from collections.abc import Sized
from dataclasses import dataclass
from typing import Generic, List, Tuple

from .base import PositionRange
from .uint_range import UIntRange
from .variant import BaseVariantT, sort_variants
from ..errors import InvalidVariantRef


@dataclass(frozen=True)
class VariantGroup(Generic[BaseVariantT], Sized):
    __slots__ = ['variants']

    variants: List[BaseVariantT]

    @classmethod
    def from_variants(cls, variants: List[BaseVariantT]):
        return cls(sort_variants(variants))

    def __len__(self) -> int:
        return len(self.variants)

    @property
    def is_empty(self):
        return len(self) == 0

    def apply(self, start: int, sequence: str, alt_offset: int = 0, ref_check: bool = False) -> Tuple[str, int]:
        alt_seq: str = sequence

        for variant in self.variants:
            offset: int = variant.start - start
            if ref_check and variant.ref_length > 0:
                exp_ref: str = sequence[offset:offset + variant.ref_length]
                if variant.ref != exp_ref:
                    raise InvalidVariantRef(
                        f"Invalid variant at {variant.genomic_position}: "
                        f"expected {variant.ref}, found {exp_ref}!")
            alt_seq = variant.mutate_from(
                alt_seq, offset + alt_offset, ref_check=False)

            # Correct following variant position
            alt_offset += variant.alt_ref_delta

        return alt_seq, alt_offset

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

    def overlaps_range(self, gr: PositionRange) -> bool:
        return any(
            r.ref_range.overlaps_range(gr, unstranded=True)
            for r in self.variants
        )

    def overlaps(self, variant: BaseVariantT) -> bool:
        return self.overlaps_range(variant.ref_range)
