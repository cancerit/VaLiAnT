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

from .base import GenomicRange, TranscriptInfo


@dataclass
class ExonInfo:
    __slots__ = ['transcript_info', 'genomic_range', 'exon_index']

    transcript_info: TranscriptInfo
    genomic_range: GenomicRange
    exon_index: int

    @property
    def transcript_id(self) -> str:
        return self.transcript_info.transcript_id

    @property
    def gene_id(self) -> str:
        return self.transcript_info.gene_id