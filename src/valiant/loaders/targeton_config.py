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

import re
from dataclasses import dataclass, replace

from .base_targeton_config import BaseTargetonConfig, parse_mutators
from ..strings.strand import Strand
from ..uint_range import UIntRange
from ..utils import get_not_none, is_in_opt_range
from .mutator_config import MutatorConfig
from .utils import parse_list, get_int_enum, parse_uint_range_from_list


CSV_HEADER = [
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'r2_start',
    'r2_end',
    'ext_vector',
    'action_vector',
    'sgrna_vector'
]


TargetonConfigField = get_int_enum('TargetonConfigField', CSV_HEADER)


# Mutation vector pattern, e.g.: `(1del), (snv, 1del), (3del)`
mutator_group_pt = r'\((\s*[\w\-_]+\s*(?:,\s*[\w\-_]+)*)?\s*\)'
mutator_vector_re = re.compile(
    r'\s*,\s*'.join([mutator_group_pt] * 3))


def parse_mutator_tuples(s: str) -> list[list[MutatorConfig]]:
    m = mutator_vector_re.match(s)

    if not m:
        raise ValueError("Invalid format for mutator vector!")

    return [
        parse_mutators(mutator_group) if mutator_group else []
        for mutator_group in m.groups()
    ]


def parse_ext_vector(s: str) -> tuple[int, int]:
    try:
        region_1_length, region_3_length = [
            int(t)
            for t in parse_list(s, n=2)
        ]
    except ValueError:
        raise ValueError("Invalid extension vector: two integers expected!")

    return region_1_length, region_3_length


@dataclass(slots=False)
class TargetonConfig(BaseTargetonConfig):
    contig: str
    strand: Strand
    region_1_length: int
    region_3_length: int
    mutators: tuple[list[MutatorConfig], list[MutatorConfig], list[MutatorConfig]]
    sgrna_ids: frozenset[str]

    def __post_init__(self) -> None:
        # Validate region 2
        if (
            (self.region_2.start < self.ref.start) or
            (self.region_2.end > self.ref.end)
        ):
            raise ValueError("Invalid region 2: exceeding targeton range!")

        # Validate region 1
        if self.region_2.start - self.region_1_length < self.ref.start:
            raise ValueError("Invalid region 1 length: exceeding targeton range!")

        # Validate region 3
        if self.region_2.end + self.region_3_length > self.ref.end:
            raise ValueError("Invalid region 3 length: exceeding targeton range!")

    def alter(self, ref_range: UIntRange, region_2_range: UIntRange) -> TargetonConfig:
        return replace(self, ref=ref_range, region_2=region_2_range)

    @property
    def base_name(self) -> str:
        return '_'.join([
            self.contig,
            str(self.ref.start),
            str(self.ref.end),
            self.strand.label
        ])

    @property
    def name(self) -> str:
        if not self.sgrna_ids:
            return self.base_name

        return '_'.join([self.base_name, *sorted(self.sgrna_ids)])

    @classmethod
    def from_list(cls, a: list[str]) -> TargetonConfig:

        # Parse mutator collections
        action_vector = a[TargetonConfigField.ACTION_VECTOR]  # type: ignore
        ma, mb, mc = parse_mutator_tuples(action_vector)

        # Parse sgRNA ID's
        sgrna_vector = a[TargetonConfigField.SGRNA_VECTOR]  # type: ignore
        sgrna_ids = frozenset(parse_list(sgrna_vector))

        def _parse_uint_range(start_field: int, end_field: int) -> UIntRange:
            return parse_uint_range_from_list(a, start_field, end_field)

        ref_range = _parse_uint_range(
            TargetonConfigField.REF_START,  # type: ignore
            TargetonConfigField.REF_END)  # type: ignore
        region_2 = _parse_uint_range(
            TargetonConfigField.R2_START,  # type: ignore
            TargetonConfigField.R2_END)  # type: ignore

        # Parse extension vector
        ext_vector = a[TargetonConfigField.EXT_VECTOR]  # type: ignore
        region_1_length, region_3_length = parse_ext_vector(ext_vector)

        return cls(
            ref_range,
            region_2,
            a[TargetonConfigField.REF_CHR],  # type: ignore
            Strand(a[TargetonConfigField.REF_STRAND]),  # type: ignore
            region_1_length,
            region_3_length,
            (ma, mb, mc),
            sgrna_ids)

    def get_region_1(self) -> UIntRange | None:
        if self.region_1_length == 0:
            return None
        return self.region_2.get_before(self.region_1_length)

    def get_region_3(self) -> UIntRange | None:
        if self.region_3_length == 0:
            return None
        return self.region_2.get_after(self.region_3_length)

    def get_const_1(self) -> UIntRange | None:
        r1 = self.get_region_1()
        if not r1 or r1.start <= self.ref.start:
            return None
        return UIntRange(self.ref.start, r1.start - 1)

    def get_const_2(self) -> UIntRange | None:
        r3 = self.get_region_3()
        if not r3 or r3.end >= self.ref.end:
            return None
        return UIntRange(r3.end + 1, self.ref.end)

    def is_in_const_1(self, pos: int) -> bool:
        return is_in_opt_range(self.get_const_1(), pos)

    def is_in_const_2(self, pos: int) -> bool:
        return is_in_opt_range(self.get_const_2(), pos)

    def get_const_regions(self) -> list[UIntRange]:
        return get_not_none([
            self.get_const_1(),
            self.get_const_2()
        ])

    def is_in_const_region(self, pos: int) -> bool:
        return self.is_in_const_1(pos) or self.is_in_const_2(pos)

    def get_region_mutators(self, i: int) -> list[MutatorConfig]:
        return self.mutators[i]

    def get_regions(self) -> list[UIntRange | None]:
        return [
            self.get_region_1(),
            self.region_2,
            self.get_region_3()
        ]

    def get_mutable_regions(self) -> list[tuple[UIntRange, list[MutatorConfig]]]:
        return [
            (r, self.get_region_mutators(i))
            for i, r in enumerate(self.get_regions())
            if r is not None and self.get_region_mutators(i)
        ]

    def get_all_regions(self) -> list[UIntRange]:
        return get_not_none([
            self.get_const_1(),
            self.get_region_1(),
            self.region_2,
            self.get_region_3(),
            self.get_const_2()
        ])
