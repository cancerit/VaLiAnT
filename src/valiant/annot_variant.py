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
from functools import lru_cache

from .codon_table import CodonTable
from .enums import MutationType
from .seq import Seq
from .strings.codon import Codon
from .strings.translation_symbol import TranslationSymbol
from .uint_range import UIntRange
from .variant import Variant


@lru_cache(maxsize=256)
def get_codon_range_offset(pos: int, frame: int = 0) -> tuple[UIntRange, int]:
    codon_offset = (pos + frame) % 3
    codon_start = pos - codon_offset
    codon_end = pos + (3 - codon_offset) - 1
    return UIntRange(codon_start, codon_end), codon_offset


@dataclass
class AnnotVariant(Variant):
    codon_offset: int
    codon_ref: Codon
    codon_alt: Codon
    aa_ref: TranslationSymbol
    aa_alt: TranslationSymbol

    def __post_init__(self) -> None:
        super().__post_init__()
        # TODO: restrict to SNV's
        if not (0 <= self.codon_offset < 3):
            raise ValueError("Invalid codon offset!")

    @classmethod
    def from_codons(cls, codon_table: CodonTable, v: Variant, offset: int, ref: Codon, alt: Codon) -> AnnotVariant:
        return cls(
            v.pos, v.ref, v.alt,
            offset, ref, alt,
            codon_table.translate(ref),
            codon_table.translate(alt))

    @property
    def codon_start(self) -> int:
        return self.pos - self.codon_offset

    def mutation_type(self) -> MutationType:
        return self.aa_ref.get_aa_change(self.aa_alt)

    @classmethod
    def annotate(cls, codon_table: CodonTable, seq: Seq, v: Variant, frame: int = 0) -> AnnotVariant:
        var_offset = seq.get_offset(v.pos)
        # TODO: handle
        assert 0 <= var_offset < len(seq)
        codon_range, codon_offset = get_codon_range_offset(
            v.pos - seq.start, frame=frame)
        # TODO: verify it is not a partial codon (expect seq to be a complete CDS?)
        codon_ref = Codon(seq.substr(codon_range, rel=True))
        codon_alt = Codon(
            codon_ref.replace_substr(UIntRange.from_pos(codon_offset), v.alt))

        return AnnotVariant.from_codons(codon_table, v, codon_offset, codon_ref, codon_alt)
