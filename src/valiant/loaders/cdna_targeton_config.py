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

from ..mutator import MutatorCollection
from ..uint_range import UIntRange
from .base_targeton_config import BaseTargetonConfig, parse_mutators
from .mutator_config import MutatorConfig
from .utils import get_int_enum, parse_uint_range_from_list


CSV_HEADER = [
    'seq_id',
    'targeton_start',
    'targeton_end',
    'r2_start',
    'r2_end',
    'action_vector'
]


TargetonConfigField = get_int_enum('TargetonConfigField', CSV_HEADER)


@dataclass(slots=False)
class CDNATargetonConfig(BaseTargetonConfig):
    seq_id: str
    mutators: frozenset[MutatorConfig]

    @property
    def name(self) -> str:
        a: list[str] = [
            self.seq_id,
            str(self.ref.start),
            str(self.ref.end)
        ]
        return '_'.join(a)

    @property
    def targeton_range(self) -> UIntRange:
        return self.ref

    @property
    def mutator_collection(self) -> MutatorCollection:
        return MutatorCollection.from_configs(list(self.mutators))

    def get_hash(self) -> str:
        from hashlib import md5

        # Used to disambiguate file names across experiments
        return md5(':'.join([
            self.seq_id,
            repr(self.targeton_range.to_tuple()),
            repr(self.region_2.to_tuple()),
            repr(sorted({
                m.as_str()
                for m in self.mutator_collection.mutators
            }))
        ]).encode('utf-8')).hexdigest()

    @classmethod
    def from_list(cls, a: list[str]) -> CDNATargetonConfig:

        # Parse sequence ID
        seq_id = a[TargetonConfigField.SEQ_ID]  # type: ignore

        def _parse_uint_range(start_field: int, end_field: int) -> UIntRange:
            return parse_uint_range_from_list(a, start_field, end_field)

        # Parse sequence ranges
        ref_range = _parse_uint_range(
            TargetonConfigField.TARGETON_START,  # type: ignore
            TargetonConfigField.TARGETON_END)  # type: ignore
        region_2 = _parse_uint_range(
            TargetonConfigField.R2_START,  # type: ignore
            TargetonConfigField.R2_END)  # type: ignore

        # Parse mutators
        mutators = frozenset(parse_mutators(
            a[TargetonConfigField.ACTION_VECTOR]))  # type: ignore

        return cls(ref_range, region_2, seq_id, mutators)
