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

from collections.abc import Sized
from dataclasses import dataclass

from .oligo_seq import alter_seq
from .seq import Seq
from .uint_range import UIntRange
from .variant import RegisteredVariant, Variant


class InvalidVariantRef(Exception):
    pass


def check_ref(ref_seq: Seq, variant: Variant) -> None:
    """Check whether the variant REF matches the reference sequence"""

    exp_ref = ref_seq.substr(variant.ref_range)
    if variant.ref != exp_ref:
        raise InvalidVariantRef(
            f"Invalid variant at {variant.pos}: "
            f"expected {variant.ref}, found {exp_ref}!")


@dataclass(frozen=True, slots=True)
class VariantGroup(Sized):
    variants: list[RegisteredVariant]

    @classmethod
    def from_variants(cls, variants: list[RegisteredVariant]):
        return cls(sorted(variants, key=lambda x: x.pos))

    def __len__(self) -> int:
        return len(self.variants)

    @property
    def is_empty(self):
        return len(self) == 0

    def apply(self, ref_seq: Seq, ref_check: bool = False) -> tuple[Seq, list[RegisteredVariant]]:
        alt_seq: Seq = ref_seq.clone()
        alt_offset: int = 0
        alt_vars: list[RegisteredVariant] = []
        offset_variant: RegisteredVariant

        for variant in self.variants:

            # Optionally validate the REF
            if ref_check and variant.ref_len > 0:
                check_ref(ref_seq, variant)

            # Offset the variant based on the preceding variants
            offset_variant = variant.offset(alt_offset)
            alt_vars.append(offset_variant)

            # Mutate the ALT sequence
            alt_seq = alter_seq(alt_seq, offset_variant)

            # Correct the following variant position
            alt_offset += variant.alt_ref_delta

        return alt_seq, alt_vars

    def get_sub(self, r: UIntRange) -> VariantGroup:
        def var_in_range(variant: Variant) -> bool:
            if variant.pos in r:
                return True
            if variant.ref_len > 1:
                var_end = variant.pos + variant.ref_len - 1
                if var_end in r:
                    return True
            return False

        return VariantGroup(list(filter(var_in_range, self.variants)))

    def overlaps_range(self, r: UIntRange) -> bool:
        return any(
            v.ref_range.intersect(r)
            for v in self.variants
        )

    def overlaps(self, variant: Variant) -> bool:
        return self.overlaps_range(variant.ref_range)
