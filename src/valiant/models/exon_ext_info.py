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
from typing import Tuple

from .base import GenomicRange, TranscriptInfo
from .exon_info import ExonInfo


@dataclass
class ExonExtInfo:
    __slots__ = [
        'exon_info',
        'strand',
        'delta_5p',
        'delta_3p',
        'cds_ext_5_length',
        'cds_ext_3_length'
    ]

    exon_info: ExonInfo
    strand: str
    delta_5p: int
    delta_3p: int
    cds_ext_5_length: int
    cds_ext_3_length: int

    @property
    def frame(self) -> int:
        return self.cds_ext_5_length

    @property
    def genomic_range(self) -> GenomicRange:
        return self.exon_info.genomic_range

    @classmethod
    def from_pyr(cls, chr_strand: Tuple[str, str], record) -> ExonExtInfo:
        return ExonExtInfo(
            ExonInfo(
                TranscriptInfo(record.gene_id, record.transcript_id),
                GenomicRange.from_pyr(chr_strand, record),
                record.exon_index
            ),
            record.Strand,
            record.delta_5p,
            record.delta_3p,
            record.frame,
            record.cds_ext_3_length
        )
