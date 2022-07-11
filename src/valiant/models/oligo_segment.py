########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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

import abc
from dataclasses import dataclass
from typing import Dict, FrozenSet, Generic, List, Optional

from valiant.enums import MutationType, TargetonMutator

from valiant.models.base import GenomicRange, StrandedPositionRange
from valiant.models.codon_table import CodonTable
from valiant.models.mutated_sequences import MutationCollection
from valiant.models.snv_table import AuxiliaryTables
from valiant.models.targeton import CDSTargeton, ITargeton, Targeton, TargetonT, PamProtected
from valiant.models.custom_variants import CustomVariant


@dataclass
class OligoSegment(abc.ABC, Generic[TargetonT]):
    region: TargetonT

    @property
    def start(self) -> int:
        return self.genomic_range.start

    @property
    def sequence(self) -> str:
        return self.region.ref_seq

    @property
    def genomic_range(self) -> GenomicRange:
        # This property can only be used for SGE libraries
        return self.region.pos_range  # type: ignore

    @property
    def pam_protected_sequence(self) -> str:
        return self.region.alt_seq

    # See mypy bug: https://github.com/python/mypy/issues/12163
    def get_pam_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        is_pam: bool = issubclass(type(self.region), PamProtected)  # type: ignore
        return (
            self.region.get_pam_variant_annotations(codon_table) if is_pam else  # type: ignore
            []
        )

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return self.region.get_sgrna_ids(spr)


@dataclass
class InvariantOligoSegment(OligoSegment, Generic[TargetonT]):
    __slots__ = ['region']


@dataclass
class TargetonOligoSegment(OligoSegment[ITargeton], Generic[TargetonT]):
    __slots__ = ['region', 'mutators']

    mutators: FrozenSet[TargetonMutator]

    @property
    def targeton(self) -> ITargeton:
        return self.region

    def compute_mutations(self, aux: AuxiliaryTables, **kwargs) -> Dict[TargetonMutator, MutationCollection]:
        if isinstance(self.targeton, CDSTargeton):
            return self.targeton.compute_mutations(self.mutators, aux, **kwargs)
        if isinstance(self.targeton, Targeton):
            return self.targeton.compute_mutations(self.mutators)
        raise TypeError("Invalid targeton type!")

    @property
    def start(self) -> int:
        return self.targeton.start
