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
from dataclasses import dataclass

from valiant.strings.strand import Strand

from .csv import load_csv
from .mutator_config import MutatorConfig
from .utils import parse_list, get_int_enum
from ..uint_range import UIntRange

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


def parse_mutators(s: str) -> list[MutatorConfig]:
    mutator_codes = sorted(set(parse_list(s)))
    return list(map(MutatorConfig.parse, mutator_codes))


def parse_mutator_tuples(s: str) -> list[list[MutatorConfig]]:
    m = mutator_vector_re.match(s)

    if not m:
        raise ValueError("Invalid format for mutator vector!")

    return [
        parse_mutators(mutator_group) if mutator_group else []
        for mutator_group in m.groups()
    ]


@dataclass(slots=True)
class TargetonConfig:
    contig: str
    strand: Strand
    ref: UIntRange
    region_2: UIntRange
    target_region_2_extension: tuple[int, int]
    mutators: tuple[list[MutatorConfig], list[MutatorConfig], list[MutatorConfig]]
    sgrna_ids: frozenset[str]

    @property
    def name(self) -> str:
        a: list[str] = [
            self.contig,
            str(self.ref.start),
            str(self.ref.end),
            self.strand.label
        ]
        if self.sgrna_ids:
            a.extend(sorted(self.sgrna_ids))
        return '_'.join(a)

    @classmethod
    def from_list(cls, a: list[str]) -> TargetonConfig:
        # Parse extension vector
        try:
            ext_vector = a[TargetonConfigField.EXT_VECTOR]  # type: ignore
            ext_a, ext_b = [
                int(t)
                for t in parse_list(ext_vector, n=2)
            ]
        except ValueError:
            raise ValueError("Invalid extension vector: two integers expected!")

        # Parse mutator collections
        action_vector = a[TargetonConfigField.ACTION_VECTOR]  # type: ignore
        ma, mb, mc = parse_mutator_tuples(action_vector)

        # Parse sgRNA ID's
        sgrna_vector = a[TargetonConfigField.SGRNA_VECTOR]  # type: ignore
        sgrna_ids = frozenset(parse_list(sgrna_vector))

        return cls(
            a[TargetonConfigField.REF_CHR],  # type: ignore
            Strand(a[TargetonConfigField.REF_STRAND]),  # type: ignore
            UIntRange(
                int(a[TargetonConfigField.REF_START]),  # type: ignore
                int(a[TargetonConfigField.REF_END])  # type: ignore
            ),
            UIntRange(
                int(a[TargetonConfigField.R2_START]),  # type: ignore
                int(a[TargetonConfigField.R2_END])  # type: ignore
            ),
            # TODO: improve validation
            (ext_a, ext_b),
            (ma, mb, mc),
            sgrna_ids)


@dataclass(slots=True)
class ExperimentConfig:
    contig: str
    strand: str
    targeton_configs: list[TargetonConfig]

    @classmethod
    def load(cls, fp: str) -> ExperimentConfig:
        targetons = [
            TargetonConfig.from_list(r)
            for r in load_csv(fp, columns=CSV_HEADER, delimiter='\t')
        ]
        if not targetons:
            raise ValueError("Invalid experiment configuration: no targetons!")

        t = targetons[0]
        return cls(t.contig, t.strand, targetons)

    @property
    def ref_ranges(self) -> list[UIntRange]:
        """Unique reference ranges"""

        return list(set(t.ref for t in self.targeton_configs))

    @property
    def sgrna_ids(self) -> frozenset[str]:
        return frozenset().union(t.sgrna_ids for t in self.targeton_configs)

    def __post_init__(self) -> None:
        if any(
            t.contig != self.contig or t.strand != self.strand
            for t in self.targeton_configs
        ):
            raise ValueError("Multiple contig and/or different strands not supported!")
