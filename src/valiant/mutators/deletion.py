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

from . import BaseMutator
from .codon import get_codon_replacements
from ..cds_seq import CdsSeq
from ..int_pattern_builder import pt_codon
from ..mutator_type import MutatorType
from ..seq import Seq
from ..variant import Variant


# TODO: add slots back when CPython fixes 90562
#  https://github.com/python/cpython/issues/90562
@dataclass(frozen=True, slots=False, init=False)
class DeletionMutator(BaseMutator):
    TYPE = MutatorType.DEL

    def as_str(self) -> str:
        # For backwards compatibility
        return (
            '1del' if self.pt.span == 1 else
            f"{self.pt.span}{self.TYPE.value}{self.pt.offset}"
        )

    def get_variants(self, seq: Seq) -> list[Variant]:
        refs = self.get_refs(seq)
        return [
            Variant.get_del(ref.start, ref.s)
            for ref in refs
        ]


@dataclass(frozen=True, slots=False, init=False)
class InFrameDeletionMutator(BaseMutator):
    TYPE = MutatorType.IN_FRAME

    def __init__(self) -> None:
        super().__init__(pt_codon)

    def get_variants(self, seq: CdsSeq) -> list[Variant]:
        return get_codon_replacements(seq, self, None)
