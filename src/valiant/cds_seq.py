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

from .seq import Seq
from .strings.codon import Codon
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .utils import get_codon_offset_complement


@dataclass(slots=False)
class CdsSeq(Seq):
    cds_prefix: DnaStr
    cds_suffix: DnaStr

    def __post_init__(self) -> None:
        assert self.ext_length % 3 == 0

    def as_codon(self) -> Codon:
        assert (
            len(self) == 3 and
            not self.cds_prefix and
            not self.cds_suffix
        )
        return Codon(self.s)

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

    def get_inner_cds_range(self) -> UIntRange:
        a = get_codon_offset_complement(self.cds_prefix_length)
        b = get_codon_offset_complement(self.cds_suffix_length)
        r = UIntRange(self.start + a, self.end - b)
        assert len(r) % 3 == 0
        return r

    def get_inner_cds_seq(self) -> Seq:
        return self.subseq(self.get_inner_cds_range(), rel=False)

    def get_codons(self) -> list[Seq]:
        return self.get_inner_cds_seq().subseq_triplets()

    def get_ext_offset(self, pos: int) -> int:
        return super().get_offset(pos) + self.cds_prefix_length

    @property
    def ext_length(self) -> int:
        return (
            self.cds_prefix_length +
            self.cds_suffix_length +
            len(self)
        )

    @classmethod
    def from_seq(cls, seq: Seq, cds_prefix: DnaStr, cds_suffix: DnaStr) -> CdsSeq:
        return cls(
            s=seq.s,
            start=seq.start,
            cds_prefix=cds_prefix,
            cds_suffix=cds_suffix)

    @property
    def ext(self) -> DnaStr:
        return DnaStr(f"{self.cds_prefix}{self.s}{self.cds_suffix}")

    def ext_substr(self, r: UIntRange, rel: bool = True) -> DnaStr:
        return self._substr(self.ext, r, rel=rel)
