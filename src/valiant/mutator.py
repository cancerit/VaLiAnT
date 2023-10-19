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
from typing import ClassVar, Type
from itertools import chain

from .mutator_type import MutatorType, parse_mutator_type
from .mutators import BaseMutator
from .mutators.deletion import DeletionMutator
from .mutators.snv import SnvMutator
from .seq import Seq
from .utils import safe_group_by
from .variant import Variant


class MutatorBuilder:
    FIXED_MUTATORS: ClassVar[dict[MutatorType, Type[BaseMutator]]] = {
        MutatorType.SNV: SnvMutator
    }

    @classmethod
    def from_type(cls, t: MutatorType) -> BaseMutator:
        assert not t.is_parametric
        try:
            return cls.FIXED_MUTATORS[t]()
        except KeyError:
            raise NotImplementedError(f"Mutator type {t} not implemented!")

    @classmethod
    def parse(cls, s: str) -> BaseMutator:
        t, params = parse_mutator_type(s)

        if t == MutatorType.DEL:
            # Parametric deletion
            assert params is not None
            offset, span = params
            return DeletionMutator(offset, span)
        else:
            # Fixed single-nucleotide or codon change
            return cls.from_type(MutatorType(s))


@dataclass(slots=True)
class MutatorCollection:
    mutators: dict[MutatorType, BaseMutator]

    @classmethod
    def from_mutators(cls, mutators: set[BaseMutator]) -> MutatorCollection:
        return cls({
            t: list(ms)
            for t, ms in safe_group_by(mutators, lambda x: x.TYPE)
        })

    def __post_init__(self) -> None:
        if (
            MutatorType.SNV_RE in self.mutators and
            MutatorType.SNV not in self.mutators
        ):
            self.mutators[MutatorType.SNV] = {
                MutatorBuilder.from_type(MutatorType.SNV)}

    def get_variants(self, seq: Seq) -> list[Variant]:
        snvs = []
        # TODO: assign a generic priority score to each variant type instead?
        if MutatorType.SNV in self.mutators:
            snvs = MutatorBuilder.from_type(MutatorType.SNV).get_variants(seq)

        variants = snvs + list(chain.from_iterable([
            m.get_variants(
                seq, snvs=snvs) if t == MutatorType.SNV_RE else m.get_variants(seq)
            for t, muts in self.mutators.items()
            for m in muts
            if t != MutatorType.SNV
        ]))

        return variants
