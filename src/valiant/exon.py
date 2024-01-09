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

from .uint_range import UIntRange
from .utils import get_cds_ext_3_length, get_codon_offset_complement


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

    def get_codon_offset(self, strand: str, pos: int) -> int:
        # TODO: verify frame convention
        return pos - self.start + (
            self.cds_prefix_length if strand == '+' else
            self.cds_suffix_length
        )

    def get_5p_3p_extensions(self, strand: str) -> tuple[int, int]:
        return (
            (self.cds_prefix_length, self.cds_suffix_length) if strand == '+' else
            (self.cds_suffix_length, self.cds_prefix_length)
        )
