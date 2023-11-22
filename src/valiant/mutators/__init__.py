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

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar

from ..annot_variant import AnnotVariant
from ..cds_seq import CdsSeq
from ..codon_table import CodonTable
from ..mutator_type import MutatorType
from ..seq import Seq
from ..uint_range import UIntRange
from ..variant import Variant
from ..int_pattern_builder import IntPatternBuilder


@dataclass(frozen=True)
class BaseMutator(ABC):
    TYPE: ClassVar[MutatorType]

    pt: IntPatternBuilder

    def as_str(self) -> str:
        return self.TYPE.value

    def get_refs(self, seq: Seq, r: UIntRange | None = None) -> list[Seq]:
        """Get the mutation start relative positions and reference sequences"""

        if r:
            if r not in seq.get_range():
                raise ValueError("Subrange out of bounds!")
            start = r.start
            length = len(r)
        else:
            start = seq.start
            length = len(seq)

        return seq.subseq_window(self.pt, start=start, length=length)

    def get_variants(self, seq: Seq) -> list[Variant]:
        raise NotImplementedError


@dataclass(frozen=True)
class BaseCdsMutator(BaseMutator, ABC):

    @abstractmethod
    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        pass

    def get_annot_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[AnnotVariant]:
        return [
            AnnotVariant.annotate(codon_table, seq, v)
            for v in self._get_variants(codon_table, seq)
        ]
