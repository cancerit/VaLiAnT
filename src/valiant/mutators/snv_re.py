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

from ..annot_variant import AnnotVariant
from ..cds_seq import CdsSeq
from ..codon_table import CodonTable
from ..constants import STOP
from ..enums import MutationType
from ..mutator_type import MutatorType
from ..strings.codon import Codon
from ..strings.dna_str import DnaStr
from ..strings.translation_symbol import TranslationSymbol
from ..variant import Variant, VariantTuple
from .codon import CodonMutator
from .snv import SnvMutator


def snv_to_snvres(codon_table: CodonTable, snv: AnnotVariant) -> list[VariantTuple]:
    def get_snvre(x: str) -> VariantTuple:
        return snv.codon_start, DnaStr(snv.codon_ref), DnaStr(x)

    def get_top_diff_codon(c: Codon, x: TranslationSymbol) -> Codon | None:
        """Return highest ranking codon other than the input codon"""

        top_codon = codon_table.get_top_codon(x)
        return (
            codon_table.get_second_best_codon(x) if top_codon == c else
            top_codon
        )

    def get_alt_codons() -> list[Codon]:
        match snv.mutation_type:
            case MutationType.SYNONYMOUS:
                return codon_table.get_synonymous_codons(snv.codon_alt)
            case MutationType.MISSENSE | MutationType.NONSENSE:
                codon = get_top_diff_codon(snv.codon_alt, snv.aa_alt)
                return [codon] if codon else []

    alts = [
        alt
        for alt in get_alt_codons()
        if alt != snv.codon_ref and alt != snv.codon_alt
    ]
    return list(map(get_snvre, alts)) if alts else []


@dataclass(frozen=True, init=False, slots=True)
class SnvReMutator(CodonMutator):
    TYPE = MutatorType.SNV_RE

    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        # TODO: avoid recomputing the SNV's
        snvs = SnvMutator().get_annot_variants(codon_table, seq)

        max_pos = seq.end - 2
        return [
            Variant.from_tuple(t)
            for t in sorted({
                snvre
                for snv in snvs
                for snvre in snv_to_snvres(codon_table, snv)
                if snvre[0] <= max_pos
            }, key=lambda t: t[0])
        ]
