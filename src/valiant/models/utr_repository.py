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
from typing import Dict

from pyranges import PyRanges

from .base import GenomicRange, TranscriptInfo


@dataclass
class UTRRepository:
    __slots__ = ['utr_ranges']

    utr_ranges: PyRanges

    def get_transcript_infos(self, ref_ranges: PyRanges) -> Dict[GenomicRange, TranscriptInfo]:
        if not ref_ranges:
            return {}

        matches: PyRanges = ref_ranges.join(self.utr_ranges, strandedness='same')

        if not matches:
            return {}

        return {
            GenomicRange.from_pyr(k, record): TranscriptInfo(
                record.gene_id,
                record.transcript_id)
            for k, ranges in matches
            for record in ranges.itertuples()
        }
