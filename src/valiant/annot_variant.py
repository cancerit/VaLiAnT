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
from .utils import replace_character_at
from .variant import Variant


@lru_cache(maxsize=256)
def get_codon_range_offset(pos: int, frame: int = 0) -> tuple[UIntRange, int]:
    codon_offset = (pos + frame) % 3
    codon_start = pos - codon_offset
    codon_end = pos + (3 - codon_offset) - 1
    return UIntRange(codon_start, codon_end), codon_offset


@dataclass
class AnnotVariant(Variant):
    codon: Codon
    codon_offset: int
    aa_ref: TranslationSymbol
    aa_alt: TranslationSymbol

    def __post_init__(self) -> None:
        super().__post_init__()
        # TODO: restrict to SNV's
        if not (0 <= self.codon_offset < 3):
            raise ValueError("Invalid codon offset!")

    @property
    def codon_start(self) -> int:
        return self.pos - self.codon_offset

    @property
    def alt_codon(self) -> Codon:
        return Codon(replace_character_at(
            self.codon, self.codon_offset, self.alt))

    @property
    def mutation_type(self) -> MutationType:
        return self.aa_ref.get_aa_change(self.aa_alt)

    @classmethod
    def annotate(cls, codon_table: CodonTable, seq: Seq, v: Variant, frame: int = 0) -> AnnotVariant:
        codon_range, codon_offset = get_codon_range_offset(v.pos, frame=frame)
        # TODO: verify it is not a partial codon (expect seq to be a complete CDS?)
        codon_seq = Codon(seq.substr(codon_range, rel=False))

        aa_ref = TranslationSymbol(codon_table.translate(codon_seq))
        aa_alt = aa_ref  # TODO: replace! DEBUG only!

        return AnnotVariant(v.pos, v.ref, v.alt, codon_seq, codon_offset, aa_ref, aa_alt)
