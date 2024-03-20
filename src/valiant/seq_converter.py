########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from typing import Sequence

from .seq import Seq
from .strings.dna_str import DnaStr
from .variant import Variant


def apply_variants(ref_seq: Seq, alt_length: int, variants: Sequence[Variant]) -> Seq:
    """
    Alter a DNA sequence based on a set of variants

    Assumptions:
    - variants are fully in range of the reference sequence
    - variants do not overlap with each other
    - variants are sorted by position
    """

    if not variants:
        return ref_seq.clone()

    # Preallocate altered sequence
    alt_seq = bytearray(alt_length)

    ref = ref_seq.s
    ref_start = ref_seq.start
    ref_length = len(ref_seq)
    n = len(variants)

    i = 0  # REF seq base offset
    j = 0  # ALT seq base offset
    k = 0  # variant index

    # Next variant (assumption: variants are sorted by position)
    v: Variant = variants[k]

    if v.pos > ref_start:

        # Copy the head of the reference sequence unaffected by variants
        delta = v.pos - ref_start
        alt_seq[:delta] = ref[:delta].encode('ascii')

        i = delta
        j = delta

    while i < ref_length:
        ref_pos = ref_start + i
        if v and v.pos == ref_pos:

            # Apply variant
            for c in v.alt:
                alt_seq[j] = ord(c)
                j += 1

            # Skip REF bases
            i += v.ref_len

            k += 1
            if k < n:

                # Step to the next variant
                v = variants[k]

            else:
                break

        else:

            # Copy REF base
            alt_seq[j] = ord(ref[i])

            i += 1
            j += 1

    if i < ref_length:

        # Copy the tail of the reference sequence unaffected by variants
        alt_seq[j:] = ref[i:].encode('ascii')

    return Seq(ref_seq.start, DnaStr(alt_seq.decode('ascii')))
