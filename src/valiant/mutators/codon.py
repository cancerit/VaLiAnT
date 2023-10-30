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

from abc import ABC, abstractmethod
from dataclasses import dataclass

from . import BaseMutator, IntPatternBuilder
from ..codon_table import CodonTable
from ..constants import STOP
from ..mutator_type import MutatorType
from ..seq import Seq
from ..strings.codon import Codon
from ..strings.dna_str import DnaStr
from ..variant import Variant

pt_codon = IntPatternBuilder(0, 3)


def get_variant_from_ref(ref: Seq, alt: DnaStr) -> Variant:
    return Variant(ref.start, ref.s, alt)


@dataclass(frozen=True, init=False)
class BaseCodonMutator(BaseMutator, ABC):
    def __init__(self) -> None:
        super().__init__(pt_codon)


@dataclass(frozen=True)
class CodonMutator(BaseCodonMutator, ABC):
    codon_table: CodonTable

    def get_top_codon(self, aa: str) -> Codon:
        return Codon(self.codon_table.get_top_codon(aa))

    def _get_codon_replacements(self, seq: Seq, value: Codon | None) -> list[Variant]:
        alt = DnaStr(value or '')
        return [
            get_variant_from_ref(ref, alt)
            for ref in self.get_refs(seq)
        ]


@dataclass(frozen=True, slots=True)
class BaseReplaceCodonMutator(CodonMutator, ABC):

    @property
    @abstractmethod
    def alt(self) -> Codon | None:
        pass

    def get_variants(self, seq: Seq) -> list[Variant]:
        return super()._get_codon_replacements(seq, self.alt)


@dataclass(frozen=True, slots=True)
class InFrameDeletionMutator(BaseReplaceCodonMutator):
    TYPE = MutatorType.IN_FRAME

    @property
    def alt(self) -> Codon | None:
        return None


@dataclass(frozen=True, slots=True)
class AlaMutator(BaseReplaceCodonMutator):
    TYPE = MutatorType.ALA

    @property
    def alt(self) -> Codon | None:
        # TODO: verify!
        return self.get_top_codon('ala')


@dataclass(frozen=True, slots=True)
class StopMutator(BaseReplaceCodonMutator):
    TYPE = MutatorType.STOP

    @property
    def alt(self) -> Codon | None:
        # TODO: verify!
        return self.get_top_codon(STOP)


@dataclass(frozen=True, slots=True)
class AminoAcidMutator(CodonMutator):
    TYPE = MutatorType.AA

    def get_alt_aa_codons(self, codon: str) -> list[str]:
        return [
            alt_codon
            for alt_codon in self.codon_table.get_codons(
                self.codon_table.translate(codon))
            if alt_codon != codon
        ]

    def get_variants(self, seq: Seq) -> list[Variant]:
        return [
            get_variant_from_ref(ref, DnaStr(alt_codon))
            for ref in self.get_refs(seq)
            for alt_codon in self.get_alt_aa_codons(ref.s)
        ]
