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
from typing import Generic, TypeVar

from .strings.dna_str import DnaStr
from .seq import Seq
from .variant import Variant


VariantT = TypeVar('VariantT', bound=Variant)


def alter_seq(seq: Seq, variant: Variant) -> Seq:
    return seq.alter(variant.ref_range, variant.alt)


@dataclass(slots=True)
class OligoSeq(Sized, Generic[VariantT]):
    variant: VariantT
    seq: DnaStr

    def __len__(self) -> int:
        return len(self.seq)

    @classmethod
    def from_ref(cls, seq: Seq, var: VariantT, rc: bool = False) -> OligoSeq:
        s = alter_seq(seq, var).s
        return cls(var, s.revc if rc else s)
