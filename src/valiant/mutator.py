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
import logging
from typing import ClassVar

from .annot_variant import AnnotVariant
from .cds_seq import CdsSeq
from .codon_table import CodonTable
from .loaders.errors import InvalidMutator
from .loaders.mutator_config import MutatorConfig
from .mutator_type import ANNOTABLE_MUTATOR_TYPES, CDS_MUTATOR_TYPES, DEPENDENT_MUTATOR_TYPES, PARAMETRIC_MUTATOR_TYPES, MutatorType
from .mutators import BaseCdsMutator, BaseMutator
from .mutators.codon import AlaMutator, StopMutator, AminoAcidMutator
from .mutators.deletion import DeletionMutator, InFrameDeletionMutator
from .mutators.snv import SnvMutator
from .mutators.snv_re import SnvReMutator
from .pattern_variant import PatternVariant
from .seq import Seq


class MutatorBuilder:

    MUTATOR_CLASSES: ClassVar[dict[MutatorType, type[BaseMutator]]] = {
        MutatorType.DEL: DeletionMutator,
        MutatorType.SNV: SnvMutator,
        MutatorType.SNV_RE: SnvReMutator,
        MutatorType.ALA: AlaMutator,
        MutatorType.STOP: StopMutator,
        MutatorType.AA: AminoAcidMutator,
        MutatorType.IN_FRAME: InFrameDeletionMutator
    }

    @classmethod
    def from_type(cls, t: MutatorType) -> BaseMutator:
        return cls.MUTATOR_CLASSES[t]()  # type: ignore

    @classmethod
    def from_config(cls, config: MutatorConfig) -> BaseMutator:
        t: MutatorType = config.type

        if t in cls.MUTATOR_CLASSES:
            mutator_cls: type[BaseMutator] = cls.MUTATOR_CLASSES[config.type]

            if config.type in PARAMETRIC_MUTATOR_TYPES:
                if not config.pt:
                    raise InvalidMutator("Parametric mutator missing parameters!")
                return mutator_cls(config.pt)

            # Assumption: all nonparametric mutator classes have zero-argument constructors
            # TODO: add abstract subclass with zero-argument constructor
            return mutator_cls()  # type: ignore

        else:
            raise NotImplementedError(f"Mutator type '{t.value}' not supported!")


@dataclass(slots=True)
class MutatorCollection:
    mutator_types: set[MutatorType]
    mutators: list[BaseMutator]

    def get_non_cds_mutators(self, exclude_annotable: bool = False) -> list[BaseMutator]:
        return [
            m for m in self.mutators
            if m.TYPE not in CDS_MUTATOR_TYPES and (
                m.TYPE not in ANNOTABLE_MUTATOR_TYPES if exclude_annotable else
                True
            )
        ]

    def get_non_annotable_cds_mutators(self) -> list[BaseMutator]:
        return [
            m for m in self.mutators
            if m.TYPE == MutatorType.IN_FRAME
        ]

    @property
    def cds_mutators(self) -> list[BaseMutator]:
        return [m for m in self.mutators if m.TYPE in CDS_MUTATOR_TYPES]

    @property
    def annotable_mutators(self) -> list[BaseCdsMutator]:
        return [
            m for m in self.mutators
            if m.TYPE in ANNOTABLE_MUTATOR_TYPES and isinstance(m, BaseCdsMutator)
        ]

    @classmethod
    def from_configs(cls, configs: list[MutatorConfig]) -> MutatorCollection:
        types: set[MutatorType] = set()
        mutators: list[BaseMutator] = []

        t: MutatorType
        for config in configs:
            t = config.type
            types.add(t)
            mutators.append(MutatorBuilder.from_config(config))

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

    def _get_variants(self, seq: Seq, mutators: list[BaseMutator]) -> list[PatternVariant]:
        return [
            PatternVariant.from_variant(m.as_str(), v)
            for m in mutators
            for v in m.get_variants(seq)
        ]

    def has_mutator_type(self, t: MutatorType) -> bool:
        return t in self.mutator_types

    def _get_annotated_variants(self, codon_table: CodonTable, seq: CdsSeq) -> list[AnnotVariant]:
        return [
            v
            for m in self.annotable_mutators
            for v in m.get_annot_variants(codon_table, seq)
        ]

    def get_variants(self, codon_table: CodonTable, seq: Seq) -> tuple[list[PatternVariant], list[AnnotVariant]]:
        if isinstance(seq, CdsSeq):
            nam = [
                *self.get_non_cds_mutators(exclude_annotable=True),
                *self.get_non_annotable_cds_mutators()
            ]
            av = self._get_annotated_variants(codon_table, seq)
        else:
            if self.cds_mutators:
                raise ValueError("CDS-only mutators assigned to non-CDS region!")
            nam = self.get_non_cds_mutators(exclude_annotable=False)
            av = []

        return self._get_variants(seq, nam), av
