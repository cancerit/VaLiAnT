########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
# 
# The Software remains the property of Genome Research Ltd ("GRL").
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
# 
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
# 
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
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
    __slots__ = {'df'}

    df: Optional[pd.DataFrame]

    def __init__(self, df: pd.DataFrame = None) -> None:
        self.df = df

    @property
    def is_empty(self) -> bool:
        return self.df is None or self.df.shape[0] == 0

    @classmethod
    def from_variants(cls, mutations):
        pass


@dataclass
class MutationCollection(BaseMutationCollection):
    __slots__ = {'df', 'mutations'}

    REQUIRED_FIELDS: ClassVar[FrozenSet[str]] = frozenset([
        'var_type',
        'mut_position',
        'ref',
        'new',
        'mseq'
    ])

    # TODO: to remove when oligonucleotide names will be computed from the metadata table
    mutations: List[MutatedSequence]

    def __init__(self, df: pd.DataFrame = None, mutations: List[MutatedSequence] = None) -> None:
        super().__init__(df)
        self.mutations = mutations or []

    def __post_init__(self) -> None:
        if self.df is not None and frozenset(self.df.columns) != self.REQUIRED_FIELDS:
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

        return cls(df=df, mutations=mutations)
