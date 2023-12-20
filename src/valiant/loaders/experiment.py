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

from ..loaders.targeton_config import TargetonConfig, CSV_HEADER
from ..strings.strand import Strand
from ..uint_range import UIntRange
from .csv import load_csv


@dataclass(slots=True)
class ExperimentConfig:
    contig: str
    strand: Strand
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
        return frozenset().union(*[t.sgrna_ids for t in self.targeton_configs])

    def __post_init__(self) -> None:
        if any(
            t.contig != self.contig or t.strand != self.strand
            for t in self.targeton_configs
        ):
            raise ValueError("Multiple contig and/or different strands not supported!")
