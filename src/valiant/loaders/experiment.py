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

from collections import defaultdict
from dataclasses import dataclass

from ..loaders.targeton_config import TargetonConfig, CSV_HEADER
from ..uint_range import UIntRange
from .csv import load_csv


def load_targetons(fp: str) -> list[TargetonConfig]:
    return [
        TargetonConfig.from_list(r)
        for r in load_csv(fp, columns=CSV_HEADER, delimiter='\t')
    ]


@dataclass(slots=True)
class ExperimentConfig:
    targeton_configs: dict[str, list[TargetonConfig]]

    @classmethod
    def from_configs(cls, targetons: list[TargetonConfig]) -> ExperimentConfig:
        if not targetons:
            raise ValueError("Invalid experiment configuration: no targetons!")

        d = defaultdict(list)
        for t in targetons:
            d[t.contig].append(t)

        return cls(d)

    @classmethod
    def load(cls, fp: str) -> ExperimentConfig:
        targetons = load_targetons(fp)
        return cls.from_configs(targetons)

    def get_contig_ref_ranges(self, contig: str) -> list[UIntRange]:
        """Unique reference ranges"""

        return list(set(
            t.ref
            for t in self.targeton_configs[contig]
        ))

    @property
    def sgrna_ids(self) -> frozenset[str]:
        return frozenset().union(*[
            t.sgrna_ids
            for tl in self.targeton_configs.values()
            for t in tl
        ])
