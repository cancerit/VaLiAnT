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

from typing import List

import numpy as np

from .variant import BaseVariantT, sort_variants


def compute_genomic_offset(variants: List[BaseVariantT]) -> int:
    """Given a list of variants, get the offset they would introduce if applied"""

    # TODO: check for overlapping variants
    return sum(
        variant.alt_ref_delta
        for variant in variants
    )


def filter_variants_by_range(ref_start: int, ref_length: int, variants: List[BaseVariantT]) -> List[BaseVariantT]:
    return [
        variant
        for variant in variants
        if 0 <= (variant.genomic_position.position - ref_start) < ref_length
    ]


def compute_alt_offsets(ref_start: int, ref_length: int, variants: List[BaseVariantT]) -> np.ndarray:
    variants_in_range: List[BaseVariantT] = sort_variants(
        filter_variants_by_range(ref_start, ref_length, variants))
    alt_length: int = ref_length + compute_genomic_offset(variants_in_range)

    alt_offsets: np.ndarray = np.zeros(alt_length, dtype=np.int32)
    alt_offset: int = 0
    for variant in variants_in_range:

        # Filter out variants that do not add to ALT
        if variant.alt_ref_delta > 0:
            ref_offset: int = variant.genomic_position.position - ref_start

            # Increment offsets spanned by ALT
            for i in range(variant.alt_length):
                alt_offset += 1
                alt_offsets[ref_offset + i] += i

            # Set all downstream offsets (assumes the variants to be sorted by position)
            alt_offsets[ref_offset + variant.alt_length:] = alt_offset

    return alt_offsets
