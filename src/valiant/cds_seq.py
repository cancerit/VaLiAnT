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

from dataclasses import dataclass, field

from .seq import Seq
from .strings.codon import Codon
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .utils import get_codon_offset_complement, is_unique_ascending


@dataclass(slots=False)
class CdsSeq(Seq):
    cds_prefix: DnaStr = DnaStr.empty()
    cds_suffix: DnaStr = DnaStr.empty()
    cds_prefix_positions: list[int] = field(default_factory=list)
    cds_suffix_positions: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert self.ext_length % 3 == 0
        assert is_unique_ascending(self.cds_prefix_positions)
        assert is_unique_ascending(self.cds_suffix_positions)

    def as_codon(self) -> Codon:
        return Codon(self.ext)

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
    def from_seq(
        cls,
        seq: Seq,
        cds_prefix: DnaStr | None = None,
        cds_suffix: DnaStr | None = None,
        cds_prefix_positions: list[int] | None = None,
        cds_suffix_positions: list[int] | None = None
    ) -> CdsSeq:
        return cls(
            s=seq.s,
            start=seq.start,
            cds_prefix=cds_prefix or DnaStr.empty(),
            cds_suffix=cds_suffix or DnaStr.empty(),
            cds_prefix_positions=cds_prefix_positions or list(),
            cds_suffix_positions=cds_suffix_positions or list())

    @property
    def ext(self) -> DnaStr:
        return DnaStr(f"{self.cds_prefix}{self.s}{self.cds_suffix}")

    def ext_substr(self, r: UIntRange, rel: bool = True) -> DnaStr:
        return self._substr(self.ext, r, rel=rel)

    @property
    def ext_positions(self) -> list[int]:
        positions = [
            *self.cds_prefix_positions,
            *self.get_range().positions,
            *self.cds_suffix_positions
        ]
        n = len(positions)
        assert n > 0 and n % 3 == 0
        # Check the positions are in ascending order and unique
        assert is_unique_ascending(positions)
        return positions

    @property
    def ext_start(self) -> int:
        # BEWARE: this may fail if the prefix positions are not set
        return self.cds_prefix_positions[0] if self.cds_prefix else self.start
