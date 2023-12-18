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

from .cds_seq import CdsSeq
from .codon_table import CodonTable
from .enums import MutationType
from .strings.codon import Codon
from .strings.translation_symbol import TranslationSymbol
from .uint_range import UIntRange
from .variant import Variant


@lru_cache(maxsize=256)
def get_codon_range_offset(pos: int) -> tuple[UIntRange, int]:
    # Assumption: the position is relative to the start of an in-frame CDS sequence
    codon_offset = pos % 3
    codon_start = pos - codon_offset
    codon_end = pos + (3 - codon_offset) - 1
    return UIntRange(codon_start, codon_end), codon_offset


@dataclass
class AnnotVariant(Variant):
    src: str
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
    def from_codons(cls, codon_table: CodonTable, v: Variant, offset: int, ref: Codon, alt: Codon, src: str = '') -> AnnotVariant:
        return cls(
            v.pos, v.ref, v.alt,
            src, offset, ref, alt,
            codon_table.translate(ref),
            codon_table.translate(alt))

    @property
    def codon_start(self) -> int:
        return self.pos - self.codon_offset

    @property
    def mutation_type(self) -> MutationType:
        return self.aa_ref.get_aa_change(self.aa_alt)

    @classmethod
    def annotate(cls, codon_table: CodonTable, seq: CdsSeq, v: Variant, src: str = '') -> AnnotVariant:
        var_offset = seq.get_ext_offset(v.pos)
        # TODO: handle
        assert 0 <= var_offset < seq.ext_length

        match v.ref_len:
            case 1:
                # SNV
                # TODO: check the frame convention
                codon_range, codon_offset = get_codon_range_offset(var_offset)
                codon_ref = seq.ext_substr(codon_range, rel=True)
                codon_alt = codon_ref.replace_substr(UIntRange.from_pos(codon_offset), v.alt)
            case 3:
                # Codon replacement
                codon_offset = 0
                codon_ref = v.ref
                codon_alt = v.alt
            case _:
                raise ValueError("Unsupported REF length for annotation!")

        return AnnotVariant.from_codons(
            codon_table,
            v,
            codon_offset,
            Codon(codon_ref),
            Codon(codon_alt),
            src=src)
