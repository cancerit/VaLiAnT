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
from collections.abc import Sized
from dataclasses import dataclass
from itertools import groupby
import logging
from typing import List, Tuple

from .base import GenomicPosition, PositionRange, GenomicRange, StrandedPositionRange
from ..errors import GenomicRangeOutOfBounds, NonUniqueError
from ..sgrna_utils import get_codon_index, get_codon_indices_in_range
from ..utils import has_duplicates, is_dna


@dataclass(frozen=True)
class Sequence(Sized):
    __slots__ = {'sequence'}

    sequence: str

    def __post_init__(self) -> None:
        if not is_dna(self.sequence):
            raise ValueError("Invalid symbols in sequence!")

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def range(self) -> PositionRange:
        return PositionRange(1, len(self))

    def get_rel_subsequence(self, pr: PositionRange) -> Sequence:
        if pr not in self.range:
            raise ValueError("Invalid subsequence: out of range!")
        return Sequence(self.sequence[pr.start - 1:pr.end])


@dataclass(frozen=True)
class ReferenceSequence(Sequence):
    __slots__ = {'sequence', 'genomic_range'}

    genomic_range: GenomicRange

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.sequence) != len(self.genomic_range):
            raise ValueError("Sequence and genomic range have different lengths!")

    def get_subsequence(self, genomic_range: GenomicRange) -> ReferenceSequence:
        start, end = self.genomic_range.get_relative_subrange(genomic_range)
        seq: str = self.sequence[start:end]
        return ReferenceSequence(seq, genomic_range)


@dataclass(frozen=True)
class CdsReferenceSequence(ReferenceSequence):
    __slots__ = {'sequence', 'genomic_range', 'cds_prefix', 'cds_suffix'}

    cds_prefix: str
    cds_suffix: str

    def __post_init__(self) -> None:
        # The sequence is guaranteed not to be an empty string
        if (
            self.cds_prefix_length > 2 or
            self.cds_suffix_length > 2 or
            self.ext_ref_seq_length % 3 != 0
        ):
            raise ValueError("Invalid extended CDS sequence length: partial codons!")

    @property
    def ext_seq_length(self) -> int:
        return self.cds_prefix_length + self.cds_suffix_length + len(self)

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

    @property
    def frame(self) -> int:
        return self.cds_prefix_length

    @property
    def codon_count(self) -> int:
        return self.ext_seq_len // 3

    @property
    def last_codon_index(self) -> int:
        return self.codon_count - 1

    @property
    def ext_seq(self) -> str:
        return f"{self.cds_prefix}{self.sequence}{self.cds_suffix}"

    def get_codon_index_at(self, p: GenomicPosition) -> int:
        return get_codon_index(self.frame, self.genomic_range.start, p)

    def get_codon_indices_at(self, ps: List[GenomicPosition], no_duplicate_codons: bool = False) -> List[int]:
        codon_indices = [self.get_ref_codon_index(x.position) for x in ps]

        # Verify no two variants affect the same codon
        if no_duplicate_codons and has_duplicates(codon_indices):
            raise NonUniqueError

        return codon_indices

    def get_codon_indices_in_range(self, pr: PositionRange) -> List[int]:
        """Get the indices of the codons spanned by the input range, if any"""

        if pr not in self.genomic_range:
            raise GenomicRangeOutOfBounds("Invalid genomic range for codon index retrieval!")

        return get_codon_indices_in_range(
            self.frame,
            self.genomic_range.start,
            self.genomic_range.end,
            pr.start,
            pr.end)
