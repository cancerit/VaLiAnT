########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from functools import partial
import logging
from typing import Dict, Iterable, Optional, Tuple
from pysam import FastaFile
from .base import GenomicRange
from ..constants import DNA_ALPHABET_SET
from ..utils import is_dna
from ..loaders.fasta import get_fasta_file


def _normalise_ref_seq(sequence: str) -> str:

    # Remove soft-masking
    norm_sequence = sequence.upper()

    # Validate sequence
    if not is_dna(norm_sequence):
        invalid_symbols = ', '.join(set(norm_sequence) - DNA_ALPHABET_SET)
        raise ValueError(f"The reference sequence contains invalid symbols: {invalid_symbols}!")

    return norm_sequence


def _fetch_sequence(fasta: FastaFile, chromosome: str, start: int, end: int) -> str:
    return fasta.fetch(reference=chromosome, start=start - 1, end=end)


@dataclass(init=False)
class ReferenceSequenceRepository:
    __slots__ = {'_sequences', '_nucleotides'}

    _sequences: Dict[str, Dict[Tuple[int, int], Optional[str]]]
    _nucleotides: Dict[str, Dict[int, Optional[str]]]

    def __init__(self) -> None:
        self._sequences = {}
        self._nucleotides = {}

    @property
    def region_count(self) -> int:
        return sum(map(len, self._sequences.values()))

    def register_region(self, chromosome: str, start: int, end: int) -> None:
        if not chromosome or start < 1 or end < 1:
            raise ValueError("Invalid genomic position!")

        t = (start, end) if end >= start else (end, start)
        if chromosome not in self._sequences:
            self._sequences[chromosome] = {t: None}
        elif t not in self._sequences[chromosome]:
            self._sequences[chromosome][t] = None

        if start > 1:
            prev_pos: int = start - 1
            if chromosome not in self._nucleotides:
                self._nucleotides[chromosome] = {prev_pos: None}
            elif prev_pos not in self._nucleotides[chromosome]:
                self._nucleotides[chromosome][prev_pos] = None

    def register_genomic_range(self, genomic_range: GenomicRange) -> None:
        self.register_region(
            genomic_range.chromosome, genomic_range.start, genomic_range.end)

    def register_genomic_ranges(self, genomic_ranges: Iterable[GenomicRange]) -> None:
        for genomic_range in genomic_ranges:
            self.register_genomic_range(genomic_range)

    def register_sequence(self, chromosome: str, start: int, end: int, sequence: str) -> None:
        # Register normalised sequence
        t = (start, end) if end >= start else (end, start)
        self._sequences[chromosome][t] = _normalise_ref_seq(sequence)

    def register_nucleotide(self, chromosome: str, position: int, nt: str) -> None:
        self._nucleotides[chromosome][position] = _normalise_ref_seq(nt)

    def fetch_sequences(self, fasta: FastaFile) -> None:
        fetch_start: int
        start_gt_1: bool
        fetch = partial(_fetch_sequence, fasta)
        logging.debug("Fetching %d reference sequences..." % self.region_count)
        for chromosome, regions in self._sequences.items():
            for start, end in regions:

                # Fetch the region sequence extended by one (unless it is at the start of the chromosome)
                start_gt_1 = start > 1
                fetch_start = (start - 1) if start > 1 else start
                logging.debug("Fetching reference sequence at %s:%d-%d." % (chromosome, fetch_start, end))
                seq: str = fetch(chromosome, fetch_start, end)

                # Register the region sequence
                self.register_sequence(chromosome, start, end, seq[1:] if start_gt_1 else seq)

                # Register nucleotide preceding the start of the region
                if start_gt_1:
                    self.register_nucleotide(chromosome, fetch_start, seq[0])

    def get_sequence(self, chromosome: str, start: int, end: int) -> Optional[str]:
        if chromosome not in self._sequences:
            return None

        t = (start, end) if end >= start else (end, start)
        return self._sequences[chromosome].get(t, None)

    def get_genomic_range_sequence(self, genomic_range: GenomicRange) -> Optional[str]:
        return self.get_sequence(
            genomic_range.chromosome, genomic_range.start, genomic_range.end)

    def get_genomic_range_subsequence(self, genomic_range: GenomicRange, start: int, end: int) -> Optional[str]:
        if start < genomic_range.start or end > genomic_range.end:
            raise ValueError("Invalid subsequence!")
        seq = self.get_genomic_range_sequence(genomic_range)
        return seq[start - genomic_range.start:end - genomic_range.start + 1] if seq else None

    def get_nucleotide(self, chromosome: str, position: int) -> Optional[str]:
        if chromosome not in self._nucleotides:
            return None

        return self._nucleotides[chromosome].get(position, None)

    def get_nucleotide_unsafe(self, chromosome: str, position: int) -> str:
        nt: Optional[str] = self.get_nucleotide(chromosome, position)
        if not nt:
            raise KeyError(f"Nucleotide at {chromosome}:{position} not available!")
        return nt


def fetch_reference_sequences(
    ref_fasta: str,
    ref_ranges: Iterable[GenomicRange]
) -> ReferenceSequenceRepository:

    # Open reference genome file
    fasta_file: FastaFile = get_fasta_file(ref_fasta)

    try:
        ref: ReferenceSequenceRepository = ReferenceSequenceRepository()
        if ref_ranges:
            ref.register_genomic_ranges(ref_ranges)
            ref.fetch_sequences(fasta_file)

    finally:

        # Close reference genome file
        fasta_file.close()

    return ref
