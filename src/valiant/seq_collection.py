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
from .uint_range import UIntRange


@dataclass
class SeqCollection:
    seqs: list[Seq]

    def __post_init__(self) -> None:
        self.seqs.sort(key=lambda s: s.start)

    def get_before(self, i: int, r: UIntRange, before: int) -> DnaStr:
        assert before >= 0

        s: DnaStr = DnaStr.empty()

        if before > 0:

            seq: Seq = self.seqs[i]
            ds = r.start - seq.start
            if before <= ds:
                # Local extension
                s += seq.substr(UIntRange(r.start - before, r.start - 1), rel=True)
            else:
                # Distal extension
                assert i > 0
                if ds > 0:
                    s += seq.head(ds)
                s += self.seqs[i - 1].tail(before - ds)

        return s

    def get_after(self, i: int, r: UIntRange, after: int) -> DnaStr:
        assert after >= 0

        s: DnaStr = DnaStr.empty()

        if after > 0:

            seq: Seq = self.seqs[i]
            ds = seq.end - r.end
            if after <= ds:
                # Local extension
                s += seq.substr(UIntRange(r.end + 1, r.end + after), rel=True)
            else:
                # Distal extension
                assert i < len(self.seqs) - 1
                s += self.seqs[i + 1].head(after)

        return s

    def get_at(self, i: int, r: UIntRange) -> DnaStr:
        return self.seqs[i].substr(r)

    def substr(self, i: int, r: UIntRange, before: int = 0, after: int = 0) -> DnaStr:
        seq: Seq = self.seqs[i]

        assert r.start >= seq.start and r.end <= seq.end

        a = self.get_before(i, r, before)
        b = self.get_at(i, r)
        c = self.get_after(i, r, after)

        return DnaStr(f"{a}{b}{c}")