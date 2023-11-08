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
from itertools import chain
import logging
from typing import ClassVar

from .loaders.errors import InvalidMutator
from .loaders.mutator_config import MutatorConfig
from .codon_table import CodonTable
from .mutator_type import DEPENDENT_MUTATOR_TYPES, PARAMETRIC_MUTATOR_TYPES, MutatorType
from .mutators import BaseMutator
from .mutators.codon import AlaMutator, StopMutator, AminoAcidMutator, InFrameDeletionMutator, CodonMutator
from .mutators.deletion import DeletionMutator
from .mutators.snv import SnvMutator
from .mutators.snv_re import SnvReMutator
from .seq import Seq
from .variant import Variant


class MutatorBuilder:

    MUTATOR_CLASSES: ClassVar[dict[MutatorType, type[BaseMutator]]] = {
        MutatorType.DEL: DeletionMutator,
        MutatorType.SNV: SnvMutator
    }

    CDS_MUTATOR_CLASSES: ClassVar[dict[MutatorType, type[CodonMutator]]] = {
        MutatorType.SNV_RE: SnvReMutator,
        MutatorType.ALA: AlaMutator,
        MutatorType.STOP: StopMutator,
        MutatorType.AA: AminoAcidMutator,
        MutatorType.IN_FRAME: InFrameDeletionMutator
    }

    @classmethod
    def _from_config_cds(cls, codon_table: CodonTable, config: MutatorConfig) -> CodonMutator:
        mutator_cls: type[CodonMutator] = cls.CDS_MUTATOR_CLASSES[config.type]
        print(mutator_cls)
        return mutator_cls(codon_table)

    @classmethod
    def _from_config_noncds(cls, config: MutatorConfig) -> BaseMutator:
        mutator_cls: type[BaseMutator] = cls.MUTATOR_CLASSES[config.type]

        if config.type in PARAMETRIC_MUTATOR_TYPES:
            if not config.pt:
                raise InvalidMutator("Parametric mutator missing parameters!")
            return mutator_cls(config.pt)

        # Assumption: all nonparametric mutator classes have zero-argument constructors
        # TODO: add abstract subclass with zero-argument constructor
        return mutator_cls()  # type: ignore

    @classmethod
    def from_type(cls, t: MutatorType) -> BaseMutator:
        return cls.MUTATOR_CLASSES[t]()  # type: ignore

    @classmethod
    def from_config(cls, codon_table: CodonTable, config: MutatorConfig) -> BaseMutator:
        t: MutatorType = config.type

        if t in cls.MUTATOR_CLASSES:
            return cls._from_config_noncds(config)

        elif t in cls.CDS_MUTATOR_CLASSES:
            return cls._from_config_cds(codon_table, config)

        else:
            raise NotImplementedError(f"Mutator type '{t.value}' not supported!")


@dataclass(slots=True)
class MutatorCollection:
    mutator_types: set[MutatorType]
    mutators: list[BaseMutator]

    @classmethod
    def from_configs(cls, codon_table: CodonTable, configs: list[MutatorConfig]) -> MutatorCollection:
        types: set[MutatorType] = set()
        mutators: list[BaseMutator] = []

        t: MutatorType
        for config in configs:
            t = config.type
            types.add(t)
            mutators.append(MutatorBuilder.from_config(codon_table, config))

        return cls(types, mutators)

    def __post_init__(self) -> None:
        for t in list(self.mutator_types):
            if t in DEPENDENT_MUTATOR_TYPES:
                for dt in DEPENDENT_MUTATOR_TYPES[t]:
                    if dt not in self.mutator_types:
                        logging.info(
                            "Adding dependent mutator '%s' (required by '%s')..." %
                            (dt.value, t.value))
                        self.mutator_types.add(dt)

                        # Assumption: no dependencies on parametric or CDS mutators
                        self.mutators.append(MutatorBuilder.from_type(dt))

    def get_variants(self, seq: Seq) -> list[Variant]:
        variants = list(chain.from_iterable([
            m.get_variants(seq)
            for m in self.mutators
        ]))

        return variants
