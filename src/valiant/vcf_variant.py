########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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

from .strings.dna_str import DnaStr
from .strings.nucleotide import Nucleotide
from .variant import Variant

first = slice(0, 1)
last = slice(-1, None)
skip_first = slice(1, None)
skip_last = slice(0, -1)


def normalise_vcf_seq(seq: DnaStr, pos_gt_one: bool = True) -> DnaStr:
    return seq.slice(skip_first if pos_gt_one else skip_last)


@dataclass(slots=True)
class VcfVariant(Variant):

    def as_variant(
        self,
        pos: int | None = None,
        ref: DnaStr | Nucleotide | None = None,
        alt: DnaStr | Nucleotide | None = None
    ) -> Variant:
        return Variant(
            pos if pos is not None else self.pos,
            DnaStr(ref) if ref is not None else self.ref,
            DnaStr(alt) if alt is not None else self.alt)

    def normalise(self) -> tuple[Variant, Nucleotide | None]:
        ref_one: bool = self.ref_len == 1
        alt_one: bool = self.alt_len == 1

        if ref_one and alt_one:

            # Substitution (SNV)
            return self.as_variant(), None

        pos_gt_one: bool = self.pos > 1

        # Collect potentially shared nucleotide
        nt_sl: slice = first if pos_gt_one else last
        ref_nt = Nucleotide(self.ref[nt_sl])
        alt_nt = Nucleotide(self.alt[nt_sl])

        if ref_nt != alt_nt:

            # Substitution (MNV)
            return self.as_variant(), None

        # Correct position
        pos = (self.pos - 1) if pos_gt_one else self.pos

        if ref_one:

            # Insertion
            ref = DnaStr.empty()
            alt = normalise_vcf_seq(self.alt, pos_gt_one=pos_gt_one)

        else:

            # Deletion
            ref = normalise_vcf_seq(self.ref, pos_gt_one=pos_gt_one)
            alt = DnaStr.empty()

        return self.as_variant(pos=pos, ref=ref, alt=alt), ref_nt
