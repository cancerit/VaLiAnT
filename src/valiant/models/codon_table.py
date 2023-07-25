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
from itertools import groupby
from typing import Callable, Dict, List, Tuple

from ..enums import MutationType
from ..globals import TRIPLET_RCS
from ..utils import validate_strand

CODON_LENGTH = 3
STOP_CODE = 'STOP'


def get_codon_mutation_type(
    translate_f: Callable[[str], str],
    ref: str,
    alt: str
) -> MutationType:
    alt_aa: str = translate_f(alt)

    if alt_aa == STOP_CODE:
        return MutationType.NONSENSE

    ref_aa: str = translate_f(ref)
    return (
        MutationType.SYNONYMOUS if alt_aa == ref_aa else
        MutationType.MISSENSE
    )


class CodonTable:
    __slots__ = {'_codon2aa', '_aa2codons', '_codon2syn', '_codonrc2syn', '_amino_acid_symbols'}

    def __init__(self, codon_table: List[Tuple[str, str, float, int]]) -> None:

        # Codon -> Amino acid
        self._codon2aa: Dict[str, str] = {
            codon: aa
            for codon, aa, _, _ in codon_table
        }

        # Amino acid -> Codons (rank sorted)
        self._aa2codons: Dict[str, List[str]] = {
            aa: [t[0] for t in sorted(ts, key=lambda x: x[3])]
            for aa, ts in groupby(codon_table, key=lambda t: t[1])
        }

        # Codon -> Synonymous codons (rank sorted)
        self._codon2syn: Dict[str, List[str]] = {
            codon: [c for c in self._aa2codons[self._codon2aa[codon]] if c != codon]
            for codon in self._codon2aa.keys()
        }

        # Codon -> Synonymous codons of reverse complement (rank sorted)
        self._codonrc2syn: Dict[str, List[str]] = {
            codon: [
                c
                for c in self._aa2codons[self._codon2aa[TRIPLET_RCS[codon]]]
                if c != TRIPLET_RCS[codon]
            ]
            for codon in self._codon2aa.keys()
        }

        self._amino_acid_symbols = sorted(self._aa2codons.keys())

    @property
    def amino_acid_symbols(self) -> List[str]:
        return self._amino_acid_symbols

    def get_top_ranking_codon_table(self) -> List[Tuple[str, str]]:
        return [(aa, codons[0]) for aa, codons in self._aa2codons.items()]

    def get_top_ranking_codon_table_rc(self) -> List[Tuple[str, str]]:
        return [(aa, TRIPLET_RCS[codons[0]]) for aa, codons in self._aa2codons.items()]

    def get_top_ranking_codon(self, aa: str) -> str:
        return self._aa2codons[aa][0]

    def get_top_ranking_codon_rc(self, aa: str) -> str:
        return TRIPLET_RCS[self._aa2codons[aa][0]]

    def translate(self, codon: str) -> str:
        return self._codon2aa[codon]

    def translate_rc(self, codon: str) -> str:
        return self._codon2aa[TRIPLET_RCS[codon]]

    def get_translate_f(self, strand: str) -> Callable[[str], str]:
        validate_strand(strand)
        return (
            self.translate if strand == '+' else
            self.translate_rc
        )

    def get_translation_table(self) -> List[Tuple[str, str]]:
        return list(self._codon2aa.items())

    def get_translation_table_rc(self) -> List[Tuple[str, str]]:
        return [
            (TRIPLET_RCS[codon], aa)
            for codon, aa in self._codon2aa.items()
        ]

    def get_mutation_types_at(
        self,
        strand: str,
        ref_seq: str,
        alt_seq: str,
        codon_indices: List[int]
    ) -> List[MutationType]:
        tr = self.get_translate_f(strand)

        def get_mutation_type(codon_index: int) -> MutationType:
            start: int = codon_index * CODON_LENGTH
            sl = slice(start, start + CODON_LENGTH)
            return get_codon_mutation_type(tr, ref_seq[sl], alt_seq[sl])

        return list(map(get_mutation_type, codon_indices))

    @classmethod
    def load(cls, fp: str) -> CodonTable:
        from ..loaders.codon_table import load_codon_table

        return cls(load_codon_table(fp))
