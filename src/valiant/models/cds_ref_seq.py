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
from typing import List

from .base import GenomicPosition, PositionRange
from .dna_str import DnaStr
from .ref_seq import RefSeq
from ..errors import GenomicRangeOutOfBounds, NonUniqueError
from ..sgrna_utils import get_codon_index, get_codon_indices_in_range
from ..utils import has_duplicates


@dataclass(frozen=True)
class CdsRefSeq(RefSeq):
    """Non-empty CDS DNA sequence with genomic coordinates"""

    __slots__ = [*RefSeq.__slots__, 'cds_prefix', 'cds_suffix']

    cds_prefix: DnaStr
    cds_suffix: DnaStr

    def __post_init__(self) -> None:
        super().__post_init__()
        # The sequence is guaranteed not to be an empty string
        if (
            self.cds_prefix_length > 2 or
            self.cds_suffix_length > 2 or
            self.ext_seq_length % 3 != 0
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
        return self.ext_seq_length // 3

    @property
    def last_codon_index(self) -> int:
        return self.codon_count - 1

    @property
    def ext_seq(self) -> str:
        return f"{self.cds_prefix}{self.sequence}{self.cds_suffix}"

    def get_codon_index_at(self, p: GenomicPosition) -> int:
        return get_codon_index(self.frame, self.genomic_range.start, p.position)

    def get_codon_indices_at(self, ps: List[GenomicPosition], no_duplicate_codons: bool = False) -> List[int]:
        codon_indices = list(map(self.get_codon_index_at, ps))

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
