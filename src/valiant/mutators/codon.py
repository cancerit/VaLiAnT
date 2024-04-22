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

from ..int_pattern_builder import pt_codon
from ..cds_seq import CdsSeq
from ..codon_table import CodonTable
from ..constants import STOP
from ..mutator_type import MutatorType
from ..seq import Seq
from ..strings.codon import Codon
from ..strings.dna_str import DnaStr
from ..strings.translation_symbol import TranslationSymbol
from ..variant import Variant
from . import BaseCdsMutator, BaseMutator
from .utils import get_variant_from_ref


# Translation symbols
ala = TranslationSymbol('A')
stop = TranslationSymbol(STOP)


def get_codon_refs(seq: CdsSeq, m: BaseMutator) -> list[Seq]:
    return m.get_refs(seq, r=seq.get_inner_cds_range())


def get_codon_replacements(seq: CdsSeq, m: BaseMutator, value: Codon | None) -> list[Variant]:
    alt = DnaStr(value or '')
    return [
        get_variant_from_ref(ref, alt)
        for ref in get_codon_refs(seq, m)
        if alt != ref.s
    ]


@dataclass(frozen=True)
class CodonMutator(BaseCdsMutator, ABC):

    def __init__(self) -> None:
        super().__init__(pt_codon)

    def _get_codon_replacements(self, seq: CdsSeq, value: Codon) -> list[Variant]:
        return get_codon_replacements(seq, self, value)


@dataclass(frozen=True, slots=False, init=False)
class BaseReplaceCodonMutator(CodonMutator, ABC):

    @abstractmethod
    def get_alt(self, codon_table: CodonTable) -> Codon:
        pass

    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        return super()._get_codon_replacements(seq, self.get_alt(codon_table))


@dataclass(frozen=True, slots=False, init=False)
class AlaMutator(BaseReplaceCodonMutator):
    TYPE = MutatorType.ALA

    def get_alt(self, codon_table: CodonTable) -> Codon:
        return codon_table.get_top_codon(ala)


@dataclass(frozen=True, slots=False, init=False)
class StopMutator(BaseReplaceCodonMutator):
    TYPE = MutatorType.STOP

    def get_alt(self, codon_table: CodonTable) -> Codon:
        return codon_table.get_top_codon(stop)


@dataclass(frozen=True, slots=False, init=False)
class AminoAcidMutator(CodonMutator):
    TYPE = MutatorType.AA

    def get_alt_aa_codons(self, codon_table: CodonTable, codon: Codon) -> list[Codon]:
        return codon_table.get_top_codons({stop, codon_table.translate(codon)})

    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        return [
            get_variant_from_ref(ref, DnaStr(alt_codon))
            for ref in get_codon_refs(seq, self)
            for alt_codon in self.get_alt_aa_codons(codon_table, Codon(ref.s))
        ]
