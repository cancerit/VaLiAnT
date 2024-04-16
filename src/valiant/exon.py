########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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

from .strings.strand import Strand
from .uint_range import UIntRange
from .utils import clamp_non_negative, get_cds_ext_3_length, get_codon_offset_complement


def get_codon_range(strand: Strand, origin: int, codon_index: int) -> UIntRange:
    codon_start_offset = 3 * codon_index
    if strand.is_plus:
        start = origin + codon_start_offset
        end = start + 2
    else:
        end = origin - codon_start_offset
        start = end - 2
    return UIntRange(start, end)


def get_codon_range_from_offset(strand: Strand, pos: int, codon_offset: int) -> UIntRange:
    if strand.is_plus:
        start = pos - codon_offset
        end = start + 2
    else:
        end = pos + codon_offset
        start = end - 2
    return UIntRange(start, end)


@dataclass(slots=True, frozen=True)
class Exon(UIntRange):
    index: int
    frame: int

    @classmethod
    def from_range(cls, r: UIntRange, index: int, frame: int) -> Exon:
        return cls(r.start, r.end, index, frame)

    @property
    def cds_prefix_length(self) -> int:
        return self.compl_frame

    @property
    def cds_suffix_length(self) -> int:
        return get_cds_ext_3_length(self.cds_prefix_length, len(self))

    @property
    def compl_frame(self) -> int:
        return get_codon_offset_complement(self.frame)

    @property
    def next_exon_frame(self) -> int:
        return self.cds_suffix_length

    @property
    def number(self) -> int:
        return self.index

    def get_first_codon_start(self, strand: Strand) -> int:
        return (
            self.start - self.cds_prefix_length if strand.is_plus else
            self.end + self.cds_prefix_length
        )

    def get_codon_index_at(self, strand: Strand, pos: int) -> int | None:
        if pos not in self:
            return None
        return clamp_non_negative(
            abs(pos - self.get_first_codon_start(strand))) // 3

    def get_codon_indices(self, strand: Strand, r: UIntRange) -> list[int]:
        t = self.intersect(r)
        if t is None:
            return []

        first = self.get_codon_index_at(strand, t.start)
        last = self.get_codon_index_at(strand, t.end)
        assert first is not None and last is not None

        return (
            list(range(first, last + 1)) if last != first else
            [first]
        )

    def get_codon_at(self, strand: Strand, pos: int) -> UIntRange | None:
        codon_index = self.get_codon_index_at(strand, pos)
        if codon_index is None:
            return None
        return self.get_codon(strand, codon_index)

    def get_codon(self, strand: Strand, codon_index: int) -> UIntRange:
        """Get the range of positions of a codon (which may be partial)"""

        assert codon_index >= 0
        origin = self.get_first_codon_start(strand)
        # Range as if the neighbouring exons were adjacent
        r = get_codon_range(strand, origin, codon_index)
        # Clamp to the range of the exon
        s = r.intersect(self)
        assert s and 1 <= len(s) <= 3
        return s

    def get_codon_offset(self, strand: Strand, pos: int) -> int:
        # TODO: verify frame convention
        return pos - self.start + (
            self.cds_prefix_length if strand.is_plus else
            self.cds_suffix_length
        )

    def get_5p_3p_extensions(self, strand: Strand) -> tuple[int, int]:
        return (
            (self.cds_prefix_length, self.cds_suffix_length) if strand.is_plus else
            (self.cds_suffix_length, self.cds_prefix_length)
        )
