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

from .base import GenomicRange
from .seq import Seq


@dataclass(frozen=True)
class RefSeq(Seq):
    """Non-empty DNA sequence with genomic coordinates"""

    __slots__ = [*Seq.__slots__, 'genomic_range']

    genomic_range: GenomicRange

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.sequence) != len(self.genomic_range):
            raise ValueError("Sequence and genomic range have different lengths!")

    def get_subsequence(self, genomic_range: GenomicRange) -> RefSeq:
        start, end = self.genomic_range.get_relative_subrange(genomic_range)
        seq = self.sequence.slice(start, end)
        return RefSeq(seq, genomic_range)