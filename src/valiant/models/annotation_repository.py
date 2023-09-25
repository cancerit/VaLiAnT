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
from typing import Optional

from pyranges import PyRanges

from .exon import CDSContextRepository
from .utr_repository import UTRRepository
from ..loaders.gff import load_gff_cds


@dataclass
class AnnotationRepository:
    __slots__ = ['cds', 'utr']

    cds: Optional[CDSContextRepository]
    utr: Optional[UTRRepository]

    @classmethod
    def from_ranges(
        cls,
        cds_ranges: Optional[PyRanges] = None,
        utr_ranges: Optional[PyRanges] = None
    ) -> AnnotationRepository:
        return cls(
            CDSContextRepository(cds_ranges) if cds_ranges is not None else None,
            UTRRepository(utr_ranges) if utr_ranges is not None else None
        )

    @classmethod
    def from_gff(cls, fp: str) -> AnnotationRepository:
        cds_ranges, utr_ranges = load_gff_cds(fp)
        return cls.from_ranges(cds_ranges, utr_ranges)
