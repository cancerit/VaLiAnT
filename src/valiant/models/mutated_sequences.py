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

from __future__ import annotations
import abc
from dataclasses import dataclass
from typing import List, Optional, Tuple, ClassVar, FrozenSet, Callable, Type
import numpy as np
import pandas as pd
from ..constants import DNA_ALPHABET_SET
from ..enums import VariantType
from ..string_mutators import delete_1, replace_single_nucleotides, delete_non_overlapping_2


@dataclass
class MutatedSequence(abc.ABC):
    __slots__ = {'position', 'sequence'}

    position: int
    sequence: str

    type: ClassVar[VariantType]

    @classmethod
    @abc.abstractmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        pass

    @property
    @abc.abstractmethod
    def ref(self) -> str:
        pass

    @property
    def new(self) -> Optional[str]:
        return None

    def to_row(self) -> Tuple[int, int, str, Optional[str], str]:
        return self.type.value, self.position, self.ref, self.new, self.sequence


@dataclass
class DeletionMutatedSequence(MutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_seq'}

    ref_seq: str

    type: ClassVar[VariantType] = VariantType.DELETION

    @property
    def ref(self) -> str:
        return self.ref_seq

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        raise NotImplementedError()


def deletion_from_sequence(
    cls: Type[DeletionMutatedSequence],
    f: Callable[[], List[Tuple[int, str, str]]]
) -> List[MutatedSequence]:
    return [cls(i, mseq, ref_seq) for i, ref_seq, mseq in f()]


@dataclass
class Deletion1MutatedSequence(DeletionMutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_seq'}

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        return [cls(i, mseq, ref_seq) for i, ref_seq, mseq in delete_1(sequence)]


@dataclass
class Deletion2Offset0MutatedSequence(DeletionMutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_seq'}

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        return deletion_from_sequence(cls, lambda: delete_non_overlapping_2(sequence, 0))


@dataclass
class Deletion2Offset1MutatedSequence(DeletionMutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_seq'}

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        return deletion_from_sequence(cls, lambda: delete_non_overlapping_2(sequence, 1))


@dataclass
class SingleNucleotideMutatedSequence(MutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_nt', 'new_nt'}

    ref_nt: str
    new_nt: str

    type: ClassVar[VariantType] = VariantType.SUBSTITUTION

    @property
    def ref(self) -> str:
        return self.ref_nt

    @property
    def new(self) -> str:
        return self.new_nt

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        return [
            cls(i, seq, ref_nt, new_nt) for i, ref_nt, new_nt, seq in replace_single_nucleotides(
                sequence, DNA_ALPHABET_SET)
        ]


@dataclass
class SingleCodonMutatedSequence(MutatedSequence):
    __slots__ = {'position', 'sequence', 'ref_codon', 'new_codon'}

    ref_codon: str
    new_codon: str

    type: ClassVar[VariantType] = VariantType.SUBSTITUTION

    @property
    def ref(self) -> str:
        return self.ref_codon

    @property
    def new(self) -> str:
        return self.new_codon

    @classmethod
    def from_sequence(cls, sequence: str) -> List[MutatedSequence]:
        raise NotImplementedError()


@dataclass
class BaseMutationCollection(abc.ABC):
    df: Optional[pd.DataFrame]

    @property
    def is_empty(self) -> bool:
        return self.df is None or self.df.shape[0] == 0

    @classmethod
    @abc.abstractmethod
    def from_variants(cls, mutations: List):
        pass


@dataclass
class MutationCollection(BaseMutationCollection):
    __slots__ = ['df']

    REQUIRED_FIELDS: ClassVar[FrozenSet[str]] = frozenset([
        'var_type',
        'mut_position',
        'ref',
        'new',
        'mseq'
    ])

    def __post_init__(self) -> None:
        if self.df is not None and not frozenset(self.df.columns) >= self.REQUIRED_FIELDS:
            raise ValueError("Invalid mutation collection fields!")

    @classmethod
    def from_variants(cls, mutations: List) -> MutationCollection:
        df: pd.DataFrame = pd.DataFrame.from_records([
            mutation.to_row()
            for mutation in mutations
        ], columns=[
            'var_type',
            'mut_position',
            'ref',
            'new',
            'mseq'
        ])

        # Compress table
        df.var_type = df.var_type.astype(np.int8)
        df.mut_position = df.mut_position.astype(np.int32)
        df.ref = df.ref.astype('category')
        df.new = df.new.astype('category')
        df.mseq = df.mseq.astype('string')

        return cls(df)
