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
from ..strings.dna_str import DnaStr
from ..strings.translation_symbol import TranslationSymbol
from ..variant import Variant
from .codon import CodonMutator
from .snv import SnvMutator


def snv_to_snvres(codon_table: CodonTable, snv: AnnotVariant) -> list[Variant]:
    def get_snvre(x: str) -> Variant:
        return Variant(snv.codon_start, DnaStr(snv.codon_ref), DnaStr(x))

    mt = snv.mutation_type
    alts = (
        codon_table.get_synonymous_codons(snv.codon_ref) if mt == MutationType.SYNONYMOUS else [
            codon_table.get_top_codon(
                snv.aa_alt if mt == MutationType.MISSENSE else
                TranslationSymbol(STOP)  # NONSENSE
            )
        ])

    return list(map(get_snvre, alts))


@dataclass(frozen=True, init=False, slots=True)
class SnvReMutator(CodonMutator):
    TYPE = MutatorType.SNV_RE

    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        # TODO: pass a MutatorType -> Variants dictionary?
        # TODO: avoid recomputing the SNV's
        snvs = SnvMutator().get_annot_variants(codon_table, seq)
        print(snvs)
        return [
            snvre
            for snv in snvs
            for snvre in snv_to_snvres(codon_table, snv)
        ]
