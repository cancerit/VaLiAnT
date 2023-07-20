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

from dataclasses import dataclass
from typing import List, Optional

from .base import GenomicRange
from .dna_str import DnaStr
from .uint_range import UIntRange
from .variant_group import VariantGroup


def is_valid_index(range_length: int, value: int) -> bool:
    return 0 <= value < range_length


@dataclass(frozen=True)
class AltSeqBuilder:
    gr: GenomicRange
    sequence: DnaStr
    variant_groups: List[VariantGroup]

    @property
    def start(self) -> int:
        return self.gr.start

    @property
    def end(self) -> int:
        return self.gr.end

    @property
    def ext_sequence(self) -> str:
        return self.sequence

    @property
    def seq_length(self) -> int:
        return len(self.sequence)

    @property
    def ext_seq_length(self) -> int:
        return self.seq_length

    def __post_init__(self) -> None:
        if self.seq_length != len(self.gr):
            raise ValueError("Mismatching position range and sequence length!")

    def get_alt(self, extend: bool = False, variant_layer: Optional[int] = None) -> str:
        n: int = len(self.variant_groups)
        last_group_index: int = (
            variant_layer if variant_layer is not None else
            n - 1
        )
        if not is_valid_index(n, last_group_index):
            raise IndexError(f"Invalid variant group index {variant_layer}!")

        alt_seq: str = self.ext_sequence if extend else self.sequence
        for g in self.variant_groups[:last_group_index + 1]:
            alt_seq = g.apply(self.start, alt_seq)

        return alt_seq

    def get_alt_length(self, extend: bool = False, variant_layer: Optional[int] = None) -> int:
        return len(self.get_alt(extend=extend, variant_layer=variant_layer))

    def get_sub(self, r: UIntRange) -> 'AltSeqBuilder':
        nr: UIntRange = r - self.start

        return AltSeqBuilder(
            r.start,
            self.sequence[nr.to_slice()],
            [
                g.get_sub(r)
                for g in self.variant_groups
            ])
