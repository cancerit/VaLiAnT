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

from .codon import CodonMutator
from .snv import SnvMutator
from ..annot_variant import AnnotVariant
from ..codon_table import CodonTable
from ..constants import STOP
from ..enums import MutationType
from ..mutator_type import MutatorType
from ..seq import Seq
from ..strings.dna_str import DnaStr
from ..variant import Variant


def snv_to_snvres(codon_table: CodonTable, snv: AnnotVariant) -> list[Variant]:
    def get_snvre(x: str) -> Variant:
        return Variant(snv.codon_start, DnaStr(snv.codon), DnaStr(x))

    mt = snv.mutation_type
    alts = (
        codon_table.get_synonymous_codons(snv.codon) if mt == MutationType.SYNONYMOUS else
        codon_table.get_top_codon(
            snv.aa_alt if mt == MutationType.MISSENSE else
            STOP  # NONSENSE
        )
    )

    return list(map(get_snvre, alts))


@dataclass(frozen=True, slots=True)
class SnvReMutator(CodonMutator):
    TYPE = MutatorType.SNV_RE

    def get_variants(self, seq: Seq) -> list[Variant]:
        # TODO: avoid recomputing the SNV's
        snvs = [
            AnnotVariant.annotate(self.codon_table, seq, snv, frame=0)
            for snv in SnvMutator().get_variants(seq)
        ]
        return [
            snvre
            for snv in snvs
            for snvre in snv_to_snvres(self.codon_table, snv)
        ]
