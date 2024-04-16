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

from dataclasses import dataclass, replace
import logging

from .cds_seq import CdsSeq
from .exon import Exon
from .genomic_position_offsets import GenomicPositionOffsets
from .seq import Seq
from .strings.strand import Strand
from .transcript_info import TranscriptInfo
from .uint_range import UIntRange, UIntRangeSortedList
from .utils import get_cds_ext_3_length


def get_range_cds_exts(strand: Strand, exon: Exon, r: UIntRange) -> tuple[int, int]:
    delta_5p = (
        (r.start - exon.start) if strand.is_plus else
        (exon.end - r.end)
    )
    assert delta_5p >= 0

    cds_prefix_length = (exon.cds_prefix_length + (delta_5p % 3)) % 3
    cds_suffix_length = get_cds_ext_3_length(cds_prefix_length, len(r))

    return (
        (cds_prefix_length, cds_suffix_length) if strand.is_plus else
        (cds_suffix_length, cds_prefix_length)
    )


def lift_exons(strand: Strand, gpo: GenomicPositionOffsets, exons: list[Exon]) -> list[Exon]:
    if len(exons) == 0:
        return []

    # Convert exon ranges
    exon_ranges = [
        new_exon
        for old_exon in exons
        if (new_exon := gpo.ref_to_alt_range(old_exon, shrink=True))
    ]

    # Sort exon ranges by exon number
    exon_ranges.sort(
        key=lambda x: x.start,
        reverse=strand.is_minus)

    # Generate the first exon, necessarily with frame zero
    new_exons = [Exon.from_range(exon_ranges[0], 0, 0)]

    prev_exon: Exon
    if len(exons) > 1:

        # Generate the other exons
        for exon_index, r in enumerate(exon_ranges[1:], start=1):
            prev_exon = new_exons[exon_index - 1]
            new_exons.append(Exon.from_range(
                r, exon_index, prev_exon.next_exon_frame))

    if len(new_exons) != len(exons):
        logging.warning(
            "The number of exons changed from %d to %d!" %
            (len(exons), len(new_exons)))

    return new_exons


@dataclass(slots=True, frozen=True)
class Transcript:
    info: TranscriptInfo
    exons: UIntRangeSortedList[Exon]

    @property
    def strand(self) -> Strand:
        return self.info.strand

    @property
    def is_empty(self) -> bool:
        return len(self.exons) == 0

    def lift_exons(self, gpo: GenomicPositionOffsets):
        return replace(self, exons=UIntRangeSortedList(lift_exons(
            self.strand, gpo, self.exons)))

    def get_exon_index(self, exon_number: int) -> int:
        if not 0 <= exon_number < len(self.exons):
            raise ValueError(f"Exon number {exon_number} not found in transcript!")

        return (
            exon_number if self.strand.is_plus else
            (len(self.exons) - exon_number - 1)
        )

    def get_exon(self, exon_number: int) -> Exon:
        exon = self.exons[self.get_exon_index(exon_number)]
        assert exon.number == exon_number
        return exon

    def get_cds_seq(self, seq: Seq, exon_number: int, r: UIntRange) -> CdsSeq:
        exon = self.get_exon(exon_number)
        return self._get_cds_seq(seq, exon, r)

    def _get_cds_seq(self, seq: Seq, exon: Exon, r: UIntRange) -> CdsSeq:
        if r not in exon:
            raise ValueError("CDS range can't span more than one exon!")

        # Compute the lengths of the CDS extensions
        before, after = get_range_cds_exts(self.strand, exon, r)
        assert 0 <= before <= 2
        assert 0 <= after <= 2

        # Fetch liminal partial exon nucleotide positions
        exon_index = self.get_exon_index(exon.number)
        before_positions = self.exons.get_before(exon_index, r, before)
        after_positions = self.exons.get_after(exon_index, r, after)

        assert len(before_positions) == before
        assert len(after_positions) == after

        # Fetch DNA sequences
        main_seq = seq.subseq(r, rel=False)
        cds_prefix = seq.get_at(before_positions)
        cds_suffix = seq.get_at(after_positions)

        return CdsSeq.from_seq(
            main_seq,
            cds_prefix=cds_prefix,
            cds_suffix=cds_suffix,
            cds_prefix_positions=before_positions,
            cds_suffix_positions=after_positions)

    def _get_codon_seq(self, seq: Seq, exon: Exon, codon_index: int) -> CdsSeq:
        codon_range = exon.get_codon(self.strand, codon_index)
        return self._get_cds_seq(seq, exon, codon_range)

    def get_codon(self, seq: Seq, exon_number: int, codon_index: int) -> CdsSeq:
        exon = self.get_exon(exon_number)
        return self._get_codon_seq(seq, exon, codon_index)

    def get_codons_in_range(self, seq: Seq, r: UIntRange) -> list[CdsSeq]:

        # Pair exons with the subranges that overlap them, if any
        exon_ranges: list[tuple[Exon, UIntRange]] = []
        for exon in self.exons:
            ri = r.intersect(exon)
            if ri:
                exon_ranges.append((exon, ri))
            elif exon_ranges:
                break

        # Collect codons
        codons = [
            self._get_codon_seq(seq, exon, codon_index)
            for exon, ri in exon_ranges
            for codon_index in exon.get_codon_indices(self.strand, ri)
        ]

        if not codons:
            return []

        # Drop duplicate codons
        unique_codons = [codons[0]]
        for i in range(1, len(codons)):
            if codons[i].ext_start != codons[i - 1].ext_start:
                unique_codons.append(codons[i])

        return unique_codons

    def get_exon_at(self, pos: int) -> Exon | None:
        for exon in self.exons:
            if pos in exon:
                return exon
        return None

    def get_codon_at(self, seq: Seq, pos: int) -> CdsSeq | None:
        exon = self.get_exon_at(pos)
        if not exon:
            return None
        codon_range = exon.get_codon_at(self.strand, pos)
        if not codon_range:
            return None
        assert pos in codon_range
        return self._get_cds_seq(seq, exon, codon_range)
