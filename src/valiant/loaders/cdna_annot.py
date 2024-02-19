########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023, 2024 Genome Research Ltd
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

from typing import Iterable

from ..exon import Exon
from ..strings.strand import Strand
from ..transcript import Transcript
from ..transcript_info import TranscriptInfo
from ..uint_range import UIntRange, UIntRangeSortedList
from .csv import load_csv, parse_str
from .utils import parse_uint_range


CSV_HEADER = [
    'seq_id',
    'gene_id',
    'transcript_id',
    'cds_start',
    'cds_end'
]


def get_faux_transcript(seq_id: str, gene_id: str | None, transcript_id: str | None, cds_range: UIntRange | None) -> Transcript:
    exons = [Exon.from_range(cds_range, 0, 0)] if cds_range else []
    return Transcript(
        TranscriptInfo(seq_id, Strand('+'), gene_id or None, transcript_id or None),
        UIntRangeSortedList(exons))


def parse_opt_uint_range(cds_start: str | None, cds_end: str | None) -> UIntRange | None:
    if bool(cds_start) != bool(cds_end):
        raise ValueError("Invalid CDS range!")
    if not cds_start or not cds_end:
        return None
    return parse_uint_range(cds_start, cds_end)


def load_annot(fp: str, ids: Iterable[str]) -> dict[str, Transcript | None]:
    return {
        seq_id: get_faux_transcript(
            parse_str(seq_id),
            parse_str(gene_id),
            parse_str(transcript_id),
            parse_opt_uint_range(cds_start, cds_end))
        for seq_id, gene_id, transcript_id, cds_start, cds_end in load_csv(
            fp, columns=CSV_HEADER, delimiter='\t')
        if seq_id in ids
    }
