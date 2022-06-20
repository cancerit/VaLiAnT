########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from typing import List, Iterable, Tuple, Optional, FrozenSet
from ..loaders.tsv import load_tsv
from ..models.base import PositionRange
from ..enums import TargetonMutator
from ..utils import parse_mutators


CSV_HEADER = [
    'seq_id',
    'targeton_start',
    'targeton_end',
    'r2_start',
    'r2_end',
    'action_vector'
]


@dataclass(frozen=True)
class CDNATargetonConfig:
    __slots__ = {'seq_id', 'targeton_range', 'r2_range', 'mutators'}

    seq_id: str
    targeton_range: PositionRange
    r2_range: PositionRange
    mutators: FrozenSet[TargetonMutator]

    def get_hash(self) -> str:
        from hashlib import md5
        return md5(':'.join([
            self.seq_id,
            repr(self.targeton_range.to_tuple()),
            repr(self.r2_range.to_tuple()),
            repr(sorted([m.value for m in self.mutators]))
        ]).encode('utf-8')).hexdigest()

    @classmethod
    def from_row(cls, row: List[str]) -> CDNATargetonConfig:
        return cls(
            row[0],
            PositionRange(int(row[1]), int(row[2])),
            PositionRange(int(row[3]), int(row[4])),
            frozenset(parse_mutators(row[5])))


@dataclass
class CDNATargetonConfigCollection:
    __slots__ = {'cts'}

    cts: List[CDNATargetonConfig]

    def __init__(self, cts: Iterable[CDNATargetonConfig]) -> None:
        self.cts = list(cts)

    def __len__(self) -> int:
        return len(self.cts)

    @property
    def mutators(self) -> FrozenSet[TargetonMutator]:
        return frozenset.union(*[
            ctc.mutators
            for ctc in self.cts
        ])

    @classmethod
    def load(cls, fp: str) -> CDNATargetonConfigCollection:
        return cls(map(CDNATargetonConfig.from_row, load_tsv(fp, CSV_HEADER)))

    @property
    def sequence_ids(self) -> FrozenSet[str]:
        return frozenset(t.seq_id for t in self.cts)
