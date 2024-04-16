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
from itertools import chain

from ..cds_seq import CdsSeq
from ..codon_table import CodonTable
from ..mutator_type import MutatorType
from ..seq import Seq
from ..strings.nucleotide import Nucleotide
from ..variant import Variant
from ..int_pattern_builder import IntPatternBuilder
from . import BaseCdsMutator

pt_nt = IntPatternBuilder(0, 1)


# TODO: add slots back when CPython fixes 90562
#  https://github.com/python/cpython/issues/90562
@dataclass(frozen=True, slots=False, init=False)
class SnvMutator(BaseCdsMutator):
    TYPE = MutatorType.SNV

    def __init__(self) -> None:
        # Would break in Python 3.11.3 if slots were set
        super().__init__(pt_nt)

    def get_variants(self, seq: Seq) -> list[Variant]:
        refs = self.get_refs(seq)
        return list(chain.from_iterable(
            Variant.get_snvs(ref.start, Nucleotide(ref.s))
            for ref in refs
        ))

    def _get_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[Variant]:
        return self.get_variants(seq)
