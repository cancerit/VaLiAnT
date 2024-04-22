########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from .cdna_targeton_config import CDNATargetonConfig, CSV_HEADER
from .csv import load_csv


@dataclass(slots=True)
class CDNAExperimentConfig:
    targeton_configs: list[CDNATargetonConfig]

    @classmethod
    def load(cls, fp: str):
        return cls([
            CDNATargetonConfig.from_list(r)
            for r in load_csv(fp, columns=CSV_HEADER, delimiter='\t')
        ])

    def __len__(self) -> int:
        return len(self.targeton_configs)

    def __post_init__(self) -> None:
        if not self.targeton_configs:
            raise ValueError("Invalid experiment configuration: no targetons!")

    @property
    def sequence_ids(self) -> frozenset[str]:
        return frozenset(t.seq_id for t in self.targeton_configs)
