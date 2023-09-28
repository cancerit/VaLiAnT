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
from typing import Tuple

from .base import GenomicRange, TranscriptInfo
from .exon_ext_info import ExonExtInfo
from .exon_info import ExonInfo


@dataclass
class CdsInfo(ExonExtInfo):
    __slots__ = [
        *ExonExtInfo.__slots__,
        'delta_5p',
        'delta_3p'
    ]

    delta_5p: int
    delta_3p: int

    @property
    def strand(self) -> str:
        return self.genomic_range.strand

    @classmethod
    def from_pyr(cls, chr_strand: Tuple[str, str], record) -> CdsInfo:
        return cls(
            TranscriptInfo(record.gene_id, record.transcript_id),
            GenomicRange.from_pyr(chr_strand, record),
            record.exon_index,
            cds_ext_5_length=record.frame,
            cds_ext_3_length=record.cds_ext_3_length,
            delta_5p=record.delta_5p,
            delta_3p=record.delta_3p
        )

    def get_exon_info(self) -> ExonInfo:
        return ExonInfo(self.transcript_info, self.genomic_range, self.exon_index)
