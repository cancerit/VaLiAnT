########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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
from .base import PositionRange, GenomicRange
from ..utils import is_dna


@dataclass(frozen=True)
class Sequence(Sized):
    __slots__ = {'sequence'}

    sequence: str

    def __post_init__(self) -> None:
        if not is_dna(self.sequence):
            raise ValueError("Invalid symbols in sequence!")

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def range(self) -> PositionRange:
        return PositionRange(1, len(self))

    def get_rel_subsequence(self, pr: PositionRange) -> Sequence:
        if pr not in self.range:
            raise ValueError("Invalid subsequence: out of range!")
        return Sequence(self.sequence[pr.start - 1:pr.end])


@dataclass(frozen=True)
class ReferenceSequence(Sequence):
    __slots__ = {'sequence', 'genomic_range'}

    genomic_range: GenomicRange

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.sequence) != len(self.genomic_range):
            raise ValueError("Sequence and genomic range have different lengths!")

    def get_subsequence(self, genomic_range: GenomicRange) -> ReferenceSequence:
        start, end = self.genomic_range.get_relative_subrange(genomic_range)
        seq: str = self.sequence[start:end]
        return ReferenceSequence(seq, genomic_range)
