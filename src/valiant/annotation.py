########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

from .uint_range import UIntRange, UIntRangeSortedList


@dataclass(slots=True, frozen=True)
class Exon(UIntRange):
    index: int
    frame: int


@dataclass(slots=True)
class Annotation:

    # Transcript locus
    contig: str
    strand: str

    # Transcript identifiers
    gene_id: str | None
    transcript_id: str | None

    # Features
    utr: UIntRangeSortedList[UIntRange]
    cds: UIntRangeSortedList[Exon]

    @property
    def cds_start(self) -> int:
        # Assumption: CDS ranges are sorted by position
        return self.cds.ranges[0].start

    @property
    def cds_end(self) -> int:
        # Assumption: CDS ranges are sorted by position
        return self.cds.ranges[-1].end

    def __post_init__(self) -> None:
        assert self.cds_start < self.cds_end
