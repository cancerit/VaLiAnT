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
from .strings.dna_str import DnaStr


@dataclass(slots=True)
class CdsSeq(Seq):
    cds_prefix: DnaStr
    cds_suffix: DnaStr

    def __post_init__(self) -> None:
        print(self)
        assert (
            self.cds_prefix_length +
            self.cds_suffix_length +
            len(self)
        ) % 3 == 0

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

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
