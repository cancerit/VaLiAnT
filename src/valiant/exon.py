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
from .utils import get_cds_ext_3_length, get_codon_offset_complement


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
    def number(self) -> int:
        return self.index

    def get_first_codon_start(self, strand: Strand) -> int:
        return (
            self.start - self.cds_prefix_length if strand.is_plus else
            self.end + self.cds_prefix_length
        )

    def get_codon_index_at(self, pos: int) -> int | None:
        return (
            ((pos - self.compl_frame) // 3) if pos in self else
            None
        )

    def get_codon_indices(self, r: UIntRange) -> list[int]:
        t = self.intersect(r)
        if t is None:
            return []

        first = self.get_codon_index_at(t.start)
        last = self.get_codon_index_at(t.end)
        assert first is not None and last is not None

        return (
            list(range(first, last + 1)) if last != first else
            [first]
        )

    def get_codon(self, strand: Strand, codon_index: int) -> UIntRange:
        origin = self.get_first_codon_start(strand)
        return get_codon_range(strand, origin, codon_index)

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
