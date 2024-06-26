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

from dataclasses import dataclass

from .codon_table_row import CodonTableRow
from .enums import MutationType
from .errors.codon_not_found import CodonNotFound
from .strings.codon import Codon
from .strings.translation_symbol import TranslationSymbol
from .utils import safe_group_by


CodonToTransl = dict[Codon, TranslationSymbol]
TranslToCodons = dict[TranslationSymbol, list[Codon]]
CodonToCodons = dict[Codon, list[Codon]]


# TODO: make it configurable for maximum generality
START_CODON = Codon('ATG')


@dataclass(slots=True, frozen=True)
class CodonTable:
    start_codon: Codon
    codon_to_aa: CodonToTransl
    aa_to_codons: TranslToCodons
    codon_to_syn: CodonToCodons

    @classmethod
    def from_list(cls, rows: list[CodonTableRow], rc: bool = False) -> CodonTable:
        start_codon = START_CODON.revc if rc else START_CODON

        rows = [r.reverse_complement() for r in rows] if rc else rows

        # Codon -> Amino acid
        codon_to_aa = {
            r.codon: r.aa
            for r in rows
        }

        # Amino acid -> Codons (rank sorted)
        aa_to_codon = {
            aa: [t.codon for t in sorted(rs, key=lambda x: x.rank)]
            for aa, rs in safe_group_by(rows, lambda x: x.aa)
        }

        # Codon -> Synonymous codons (rank sorted)
        codon_to_syn = {
            codon: [c for c in aa_to_codon[codon_to_aa[codon]] if c != codon]
            for codon in codon_to_aa.keys()
        }

        return cls(start_codon, codon_to_aa, aa_to_codon, codon_to_syn)

    def get_codons(self, aa: TranslationSymbol) -> list[Codon]:
        return self.aa_to_codons[aa]

    def get_synonymous_codons(self, codon: Codon) -> list[Codon]:
        aa = self.translate(codon)
        return [
            syn_codon
            for syn_codon in self.get_codons(aa)
            if syn_codon != codon
        ]

    def get_codon_at(self, aa: TranslationSymbol, i: int) -> Codon | None:
        assert i >= 0
        codons = self.get_codons(aa)
        return codons[i] if i < len(codons) else None

    def get_second_best_codon(self, aa: TranslationSymbol) -> Codon | None:
        """Get second highest-ranking (most frequent) codon translating to an amino acid"""

        # Assumption: aa_to_codons has its values (codons) sorted by rank
        return self.get_codon_at(aa, 1)

    def get_top_codon(self, aa: TranslationSymbol) -> Codon:
        """Get highest-ranking (most frequent) codon translating to an amino acid"""

        # Assumption: aa_to_codons has its values (codons) sorted by rank
        return self.get_codons(aa)[0]

    def get_top_codons(self, exclude: set[TranslationSymbol] | None = None) -> list[Codon]:
        exclude = exclude or set()
        return sorted([
            self.get_top_codon(aa)
            for aa in self.aa_to_codons.keys()
            if aa not in exclude
        ])

    def translate(self, codon: Codon) -> TranslationSymbol:
        try:
            return self.codon_to_aa[codon]
        except KeyError:
            raise CodonNotFound(codon)

    def get_aa_change(self, codon_a: Codon, codon_b: Codon) -> MutationType:
        aa_a = self.translate(codon_a)
        aa_b = self.translate(codon_b)

        # Compare the amino acids
        return aa_a.get_aa_change(aa_b)

    def is_syn(self, ref: Codon, alt: Codon) -> bool:
        return self.get_aa_change(ref, alt) == MutationType.SYNONYMOUS
