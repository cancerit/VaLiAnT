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
from enum import IntEnum, Enum
from itertools import chain

from .csv import load_csv, parse_opt_str
from ..annotation import Annotation
from ..exon import Exon
from ..strings.strand import Strand
from ..transcript_info import TranscriptInfo
from ..uint_range import UIntRange, UIntRangeSortedList


class GtfField(IntEnum):
    SEQNAME = 0
    SOURCE = 1
    FEATURE = 2
    START = 3
    END = 4
    SCORE = 5
    STRAND = 6
    FRAME = 7
    ATTRIBUTE = 8


class GtfAttribute(Enum):
    GENE_ID = 'gene_id'
    TRANSCRIPT_ID = 'transcript_id'
    EXON_NUMBER = 'exon_number'


class GtfFeatureType(Enum):
    CDS = 'CDS'
    UTR = 'UTR'


def parse_gtf_attributes(s: str) -> dict[str, str]:
    return dict(
        a.strip().split(' ', 1)
        for a in s.split(';') if a
    ) if s else {}


@dataclass(slots=True, frozen=True)
class GtfFeature(UIntRange):
    gene_id: str | None
    transcript_id: str | None

    @classmethod
    def from_gtf(cls, r: list[str]) -> GtfFeature:
        attrs = parse_gtf_attributes(r[GtfField.ATTRIBUTE])
        return cls(
            int(r[GtfField.START]),
            int(r[GtfField.END]),
            parse_opt_str(attrs.get(GtfAttribute.GENE_ID.value)),
            parse_opt_str(attrs.get(GtfAttribute.TRANSCRIPT_ID.value)))

    def to_range(self) -> UIntRange:
        return UIntRange(self.start, self.end)


UtrFeature = GtfFeature


@dataclass(slots=True, frozen=True)
class CdsFeature(GtfFeature):
    frame: int

    def __post_init__(self) -> None:
        if not (0 <= self.frame <= 2):
            raise ValueError("Invalid frame in GTF/GFF2 file: out of range!")

    @classmethod
    def from_gtf(cls, r: list[str]) -> CdsFeature:
        try:
            frame = int(r[GtfField.FRAME])
        except ValueError:
            raise ValueError("Invalid frame in GTF/GFF2 file: not an integer!")

        # TODO: consider using exon_number as well? Is it guaranteed to be available?
        f = GtfFeature.from_gtf(r)
        return cls(f.start, f.end, f.gene_id, f.transcript_id, frame)

    def to_exon(self, index: int) -> Exon:
        return Exon(self.start, self.end, index, self.frame)


def cds_features_to_exons(strand: Strand, cds: list[CdsFeature]) -> UIntRangeSortedList[Exon]:

    # Order by prospective exon number
    fts = sorted(cds, key=lambda x: x.start, reverse=strand.is_minus)

    # Add stop codon to the last exon
    last_exon_index = -1
    fts[last_exon_index] = (
        fts[last_exon_index].offset_end(3) if strand.is_plus else
        fts[last_exon_index].offset_start(-3)
    )

    # Assign exon numbers
    return UIntRangeSortedList([
        ft.to_exon(i)
        for i, ft in enumerate(fts)
    ])


@dataclass(slots=True)
class GtfLoader:
    contig: str
    strand: Strand

    def load_gtf(self, fp: str) -> Annotation:
        utr: list[GtfFeature] = []
        cds: list[CdsFeature] = []

        for r in load_csv(fp, delimiter='\t'):
            # Filter by contig and strand
            if (
                r[GtfField.SEQNAME] == self.contig and
                r[GtfField.STRAND] == self.strand
            ):
                # Filter by type of feature
                match r[GtfField.FEATURE]:
                    case GtfFeatureType.CDS.value:
                        cds.append(CdsFeature.from_gtf(r))
                    case GtfFeatureType.UTR.value:
                        utr.append(UtrFeature.from_gtf(r))

        # Validate transcript identifiers
        id_pair = {
            (f.gene_id, f.transcript_id)
            for f in chain.from_iterable([utr, cds])
        }
        if len(id_pair) > 1:
            raise ValueError("Multiple transcript in annotation: not supported!")

        gene_id, transcript_id = list(id_pair)[0]
        utr_ranges = UIntRangeSortedList([x.to_range() for x in utr])

        # Assign exon numbers and add final stop codon
        exons = cds_features_to_exons(self.strand, cds)
        exons.sort()

        return Annotation(
            TranscriptInfo(
                self.contig,
                self.strand,
                gene_id,
                transcript_id),
            utr_ranges, exons)
