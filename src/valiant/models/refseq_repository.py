########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
# 
# The Software remains the property of Genome Research Ltd ("GRL").
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
# 
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
# 
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
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


# TODO: use explicit index path instead?
def get_fasta_file(fp: str) -> FastaFile:
    try:
        return FastaFile(fp)
    except IOError as ex:
        logging.critical("Failed to load reference file!")
        raise ex


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
