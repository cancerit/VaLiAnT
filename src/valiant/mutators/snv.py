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

from . import SingleBaseMutator
from ..mutator_type import MutatorType
from ..seq import Seq
from ..strings.nucleotide import Nucleotide
from ..variant import Variant


@dataclass(frozen=True, slots=True)
class SnvMutator(SingleBaseMutator):
    TYPE = MutatorType.SNV

    def get_variants(self, seq: Seq) -> list[Variant]:
        refs = self.get_refs(seq)
        return list(chain.from_iterable(
            Variant.get_snvs(ref.start, Nucleotide(ref.s))
            for ref in refs
        ))
