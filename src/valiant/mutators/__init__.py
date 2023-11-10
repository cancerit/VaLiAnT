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

from ..mutator_type import MutatorType
from ..seq import Seq
from ..variant import Variant
from .int_pattern_builder import IntPatternBuilder


@dataclass(frozen=True)
class BaseMutator(ABC):
    TYPE: ClassVar[MutatorType]

    pt: IntPatternBuilder

    def as_str(self) -> str:
        return self.TYPE.value

    def get_refs(self, seq: Seq) -> list[Seq]:
        """Get the mutation start relative positions and reference sequences"""

        starts = self.pt.build(seq.start, len(seq) - 1)
        return [
            seq.subseq(self.pt.get_range(start), rel=False)
            for start in starts
        ]

    @abstractmethod
    def get_variants(self, seq: Seq) -> list[Variant]:
        pass
