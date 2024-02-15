########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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
from typing import Generator, Generic, Iterable

from .oligo_seq import alter_seq
from .seq import Seq
from .uint_range import UIntRange
from .variant import Variant, VariantT, compute_genomic_offset, variant_sort_key


class InvalidVariantRef(Exception):
    pass


def check_ref(ref_seq: Seq, variant: Variant) -> None:
    """Check whether the variant REF matches the reference sequence"""

    exp_ref = ref_seq.substr(variant.ref_range, rel=False)
    if variant.ref != exp_ref:
        raise InvalidVariantRef(
            f"Invalid variant at {variant.pos}: "
            f"expected {variant.ref}, found {exp_ref}!")


@dataclass(slots=True, init=False)
class VariantGroup(Sized, Generic[VariantT]):
    variants: list[VariantT]

    def __init__(self, variants: Iterable[VariantT]) -> None:
        self.variants = sorted(variants, key=variant_sort_key)

    def __len__(self) -> int:
        return len(self.variants)

    @property
    def is_empty(self):
        return len(self) == 0

    def apply_no_offset(self, ref_seq: Seq) -> Seq:
        alt_seq: Seq = ref_seq.clone()
        for variant in self.variants:
            alt_seq = alter_seq(alt_seq, variant)
        return alt_seq

    def iter_alt_variants(self, alt_offset: int | None = None) -> Generator[VariantT, None, None]:
        alt_offset = alt_offset or 0
        for variant in self.variants:
            yield variant.offset(alt_offset)
            alt_offset += variant.alt_ref_delta

    def get_alt_length(self, ref_length: int) -> int:
        return ref_length + compute_genomic_offset(self.variants)

    def apply(self, ref_seq: Seq, start: int | None = None, ref_check: bool = False) -> Seq:
        alt_seq = ref_seq.clone()
        alt_offset = 0

        for variant in self.variants:

            # Optionally validate the REF
            if ref_check and variant.ref_len > 0:
                check_ref(ref_seq, variant)

            # Offset the variant based on the preceding variants
            offset_variant = variant.offset(alt_offset)

            # Mutate the ALT sequence
            alt_seq = alter_seq(alt_seq, offset_variant)

            # Correct the following variant position
            alt_offset += variant.alt_ref_delta

        if start is not None:

            # Override the start position
            alt_seq.start = start

        return alt_seq

    def overlaps_range(self, r: UIntRange) -> bool:
        return any(
            v.ref_range.intersect(r)
            for v in self.variants
        )

    def overlaps(self, variant: Variant) -> bool:
        return self.overlaps_range(variant.ref_range)
