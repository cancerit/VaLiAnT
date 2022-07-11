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
from dataclasses import dataclass, field
from functools import partial
from typing import Any, Callable, ClassVar, Dict, Generic, List, FrozenSet, Optional, Type, TypeVar
import numpy as np
import pandas as pd

from valiant.models.annotated_sequence import AnnotatedSequencePair, AnnotatedSequenceT, CDSAnnotatedSequencePair, RangeT, VariantT
from valiant.models.base import GenomicRange, StrandedPositionRange
from valiant.models.pam_protection import PamProtectedReferenceSequence, PamVariant
from .codon_table import CodonTable, STOP_CODE
from ..constants import GENERIC_MUTATORS, CDS_ONLY_MUTATORS, META_MUT_POSITION, META_PAM_MUT_SGRNA_ID, META_REF, ARRAY_SEPARATOR
from .mutated_sequences import (
    DeletionMutatedSequence,
    Deletion1MutatedSequence,
    Deletion2Offset0MutatedSequence,
    Deletion2Offset1MutatedSequence,
    MutationCollection,
    SingleNucleotideMutatedSequence,
    SingleCodonMutatedSequence
)
from .snv_table import AuxiliaryTables
from ..enums import MutationType, TargetonMutator, VariantType
from ..sgrna_utils import set_metadata_sgrna_ids
from ..string_mutators import delete_non_overlapping_3_offset, replace_codons_const
from ..utils import get_constant_category, get_out_of_frame_offset


TargetonT = TypeVar('TargetonT', bound='ITargeton')
NonCDSTargetonT = TypeVar('NonCDSTargetonT', bound='Targeton')
CDSTargetonT = TypeVar('CDSTargetonT', bound='CDSTargeton')


def get_snv_mutations(sequence: str) -> MutationCollection:
    return MutationCollection.from_variants(
        SingleNucleotideMutatedSequence.from_sequence(
            sequence))


class BaseTargeton(abc.ABC):
    MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = GENERIC_MUTATORS

    @property
    @abc.abstractmethod
    def sequence(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def pam_seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def start(self) -> int:
        pass

    def _get_mutator_method(self, mutator: TargetonMutator):
        return getattr(self, f"get_{mutator.value.replace('-', '_')}_mutations")

    def get_mutations(self, mutator: TargetonMutator, aux_tables: AuxiliaryTables = None) -> MutationCollection:
        if mutator not in self.MUTATORS:
            raise ValueError(f"Invalid mutator '{mutator.value}' for targeton!")

        return self._get_mutator_method(mutator)(aux_tables=aux_tables)

    def get_1del_mutations(self, **kwargs) -> MutationCollection:
        return MutationCollection.from_variants(
            Deletion1MutatedSequence.from_sequence(self.sequence))

    def get_2del0_mutations(self, **kwargs) -> MutationCollection:
        return MutationCollection.from_variants(
            Deletion2Offset0MutatedSequence.from_sequence(self.sequence))

    def get_2del1_mutations(self, **kwargs) -> MutationCollection:
        return MutationCollection.from_variants(
            Deletion2Offset1MutatedSequence.from_sequence(self.sequence))

    def get_snv_mutations(self, **kwargs) -> MutationCollection:
        return get_snv_mutations(self.sequence)

    def _compute_mutations(
        self,
        mutators: FrozenSet[TargetonMutator],
        aux_tables: AuxiliaryTables = None
    ) -> Dict[TargetonMutator, MutationCollection]:
        return {
            mutator: self.get_mutations(mutator, aux_tables=aux_tables)
            for mutator in mutators
        }


@dataclass(frozen=True)
class ITargeton(BaseTargeton, Generic[AnnotatedSequenceT], abc.ABC):
    annotated_seq: AnnotatedSequenceT

    @property
    def ref_seq(self) -> str:
        return self.annotated_seq.ref_seq

    @property
    def alt_seq(self) -> str:
        return self.annotated_seq.alt_seq

    @property
    def pos_range(self) -> StrandedPositionRange:
        return self.annotated_seq.pos_range

    @property
    def variant_count(self) -> int:
        return self.annotated_seq.variant_count

    @abc.abstractmethod
    def get_codon_indices(self, spr: StrandedPositionRange) -> List[int]:
        pass

    @abc.abstractmethod
    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        pass


@dataclass(frozen=True)
class Targeton(ITargeton[AnnotatedSequencePair], Generic[VariantT, RangeT]):
    __slots__ = ['annotated_seq']

    @property
    def sequence(self) -> str:
        return self.alt_seq

    @property
    def seq(self) -> str:
        return self.ref_seq

    @property
    def pam_seq(self) -> str:
        return self.alt_seq

    @property
    def start(self) -> int:
        return self.pos_range.start

    @classmethod
    def build_without_variants(
        cls: Type[NonCDSTargetonT],
        pos_range: RangeT,
        ref_seq: str
    ) -> NonCDSTargetonT:
        return cls(AnnotatedSequencePair(pos_range, ref_seq, ref_seq, []))

    @classmethod
    def from_pam_seq(
        cls: Type[NonCDSTargetonT],
        pam_seq: PamProtectedReferenceSequence
    ) -> NonCDSTargetonT:
        return cls(AnnotatedSequencePair(
            pam_seq.genomic_range,
            pam_seq.sequence,
            pam_seq.pam_protected_sequence,
            pam_seq.pam_variants))

    def get_codon_indices(self, spr: StrandedPositionRange) -> List[int]:
        return []

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return frozenset()

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator]) -> Dict[TargetonMutator, MutationCollection]:
        return super()._compute_mutations(mutators)


@dataclass(frozen=True)
class CDSTargeton(ITargeton[CDSAnnotatedSequencePair], Generic[VariantT, RangeT]):
    __slots__ = ['annotated_seq']

    MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = GENERIC_MUTATORS | CDS_ONLY_MUTATORS

    SNVRE_MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = frozenset([
        TargetonMutator.SNV_RE
    ])

    @classmethod
    def build_without_variants(
        cls: Type[CDSTargetonT],
        pos_range: RangeT,
        ref_seq: str,
        cds_prefix: str,
        cds_suffix: str
    ) -> CDSTargetonT:
        return cls(CDSAnnotatedSequencePair(
            pos_range, ref_seq, ref_seq, [], cds_prefix, cds_suffix))

    @classmethod
    def from_pam_seq(
        cls: Type[CDSTargetonT],
        pam_seq: PamProtectedReferenceSequence,
        cds_prefix: str,
        cds_suffix: str
    ) -> CDSTargetonT:
        return cls(CDSAnnotatedSequencePair(
            pam_seq.genomic_range,
            pam_seq.sequence,
            pam_seq.pam_protected_sequence,
            pam_seq.pam_variants,
            cds_prefix,
            cds_suffix))

    @property
    def sequence(self) -> str:
        return self.alt_seq

    @property
    def seq(self) -> str:
        return self.ref_seq

    @property
    def pam_seq(self) -> str:
        return self.alt_seq

    @property
    def strand(self) -> str:
        return self.pos_range.strand

    @property
    def cds_sequence(self) -> str:
        return self.annotated_seq.ext_alt_seq

    @property
    def cds_prefix_length(self) -> int:
        return self.annotated_seq.cds_prefix_length

    @property
    def cds_suffix_length(self) -> int:
        return self.annotated_seq.cds_suffix_length

    @property
    def start(self) -> int:
        return self.pos_range.start

    @property
    def frame(self) -> int:
        return self.annotated_seq.frame

    @property
    def cds_sequence_start(self) -> int:
        return self.start - self.frame

    def get_codon_indices(self, spr: StrandedPositionRange) -> List[int]:
        """Get the indices of the codons spanned by the input range, if any"""

        # Check an intersection exists
        if spr not in self.annotated_seq.pos_range:
            # Fail on incorrect strand (pathological state)
            if spr.strand != self.annotated_seq.pos_range.strand:
                raise ValueError("Failed to retrieve codon indices: incorrect strand!")
            return []

        start: int = self.annotated_seq.pos_range.start
        end: int = self.annotated_seq.pos_range.end

        # Get first and last codon indices
        codon_start: int = 0 if spr.start <= start else (spr.start - start) // 3
        codon_end: int = (
            self.annotated_seq.last_codon_index if spr.end >= end else
            self.annotated_seq.last_codon_index - ((end - spr.end) // 3)
        )

        # Generate full codon index range
        return list(range(codon_start, codon_end + 1))

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return frozenset()

    def _add_snv_metadata(
        self,
        snv_meta_full: pd.DataFrame,
        snvs: MutationCollection
    ) -> MutationCollection:

        if snvs.df is None:
            return snvs

        # Make SNV data comparable with pre-computed metadata
        snv_meta = snv_meta_full[[
            'pos',
            'alt',
            'ref_aa',
            'alt_aa',
            'mut_type'
        ]].rename(columns={
            'pos': 'mut_position',
            'alt': 'new',
            'ref_aa': 'ref_aa',
            'alt_aa': 'alt_aa',
            'mut_type': 'mut_type'
        }).set_index(['mut_position', 'new'])
        snvs.df.mut_position += self.start
        df = snvs.df.set_index(['mut_position', 'new'], drop=False)

        # Join SNV data with pre-computed metadata
        snv_joint = df.join(snv_meta)
        snv_joint = snv_joint.reset_index(drop=True)

        # TODO: to remove once the later step do not update the position
        snv_joint.mut_position -= self.start

        # Wrap complete SNV metadata in a collection
        return MutationCollection(snv_joint)

    def _get_snvres(self, aux: AuxiliaryTables, snvs: pd.DataFrame) -> MutationCollection:
        df: pd.DataFrame = aux.snvre_table.get_snvres(
            self.pos_range, self.frame, self.sequence, snvs).rename(
                columns={
                    'pos': 'mut_position',
                    'alt': 'new',
                    'ref_aa': 'ref_aa',
                    'alt_aa': 'alt_aa',
                    'mut_type': 'mut_type'
                })

        return MutationCollection(df)

    def _get_codon_mutations(self, codon_table: CodonTable, aa: str) -> MutationCollection:

        # Retrieve top codon for amino acid (or stop)
        get_top_ranking_codon: Callable[[str], str] = (
            codon_table.get_top_ranking_codon if self.strand == '+' else
            codon_table.get_top_ranking_codon_rc
        )
        codon: str = get_top_ranking_codon(aa)

        # Generate single-codon substitutions
        offset: int = get_out_of_frame_offset(self.frame)
        mc: MutationCollection = MutationCollection.from_variants([
            SingleCodonMutatedSequence(pos, mseq, ref, alt)
            for pos, ref, alt, mseq in replace_codons_const(self.sequence, offset, codon)
        ])

        if mc.df is None:
            raise RuntimeError("Mutation collection table not found!")

        # Get amino acid information
        tr: Callable[[str], str] = (
            codon_table.translate if self.strand == '+' else
            codon_table.translate_rc
        )
        rown: int = mc.df.shape[0]
        amino_acid_symbols: List[str] = codon_table.amino_acid_symbols
        mc.df['ref_aa'] = pd.Categorical(mc.df.ref.apply(tr), categories=amino_acid_symbols)
        mc.df['alt_aa'] = get_constant_category(aa, rown, amino_acid_symbols)

        # Add variant type
        mc.df['var_type'] = np.int8(VariantType.SUBSTITUTION.value)

        return mc

    def get_ala_mutations(self, aux_tables: AuxiliaryTables = None) -> MutationCollection:
        if not aux_tables:
            raise RuntimeError("Codon table not provided!")
        codon_table: CodonTable = aux_tables.codon_table
        return self._get_codon_mutations(codon_table, 'A')

    def get_stop_mutations(self, aux_tables: AuxiliaryTables = None) -> MutationCollection:
        if not aux_tables:
            raise RuntimeError("Codon table not provided!")
        codon_table: CodonTable = aux_tables.codon_table
        return self._get_codon_mutations(codon_table, STOP_CODE)

    def get_aa_mutations(self, aux_tables: AuxiliaryTables = None) -> MutationCollection:
        if not aux_tables:
            raise RuntimeError("Auxiliary tables not provided!")
        df: pd.DataFrame = aux_tables.all_aa_table.get_subs(
            self.pos_range, self.frame, self.sequence)
        return MutationCollection(df)

    def compute_mutations(
        self,
        mutators: FrozenSet[TargetonMutator],
        aux: AuxiliaryTables,
        **kwargs
    ) -> Dict[TargetonMutator, MutationCollection]:

        # Classify mutators
        base_mutators = mutators - self.SNVRE_MUTATORS
        snvre_mutators = mutators & self.SNVRE_MUTATORS

        if snvre_mutators:
            base_mutators |= {TargetonMutator.SNV}

        # Compute base mutations (overwrites existing)
        mutations: Dict[TargetonMutator, MutationCollection] = super()._compute_mutations(
            base_mutators, aux_tables=aux)

        if TargetonMutator.SNV in base_mutators:

            # Retrieve pre-computed metadata
            snv_meta_full = aux.snv_table.get_snvs(
                self.strand,
                self.cds_sequence,
                self.pos_range,
                self.cds_prefix_length,
                self.cds_suffix_length,
                reset_index=False)

            # Attach pre-computed metadata to SNV's
            snvs: MutationCollection = self._add_snv_metadata(snv_meta_full, mutations[TargetonMutator.SNV])
            mutations[TargetonMutator.SNV] = snvs

            # TODO: refactor dispatch
            for mutator in snvre_mutators:
                if mutator == TargetonMutator.SNV_RE:

                    # Compute SNVRE mutations
                    mutations[TargetonMutator.SNV_RE] = self._get_snvres(aux, snv_meta_full)

                else:
                    raise NotImplementedError(f"Mutator '{mutator.value}' not yet implemented!")

        return mutations

    def get_inframe_mutations(self, **kwargs) -> MutationCollection:
        start_offset: int = get_out_of_frame_offset(self.cds_prefix_length)
        end_offset: int = get_out_of_frame_offset(self.cds_suffix_length)
        return MutationCollection.from_variants([
            DeletionMutatedSequence(pos, mseq, ref_seq)
            for pos, ref_seq, mseq in delete_non_overlapping_3_offset(
                self.sequence, start_offset, end_offset)
        ])


class PamProtected(abc.ABC):

    @abc.abstractmethod
    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        pass


@dataclass(frozen=True)
class PamProtTargeton(Targeton[PamVariant, GenomicRange], PamProtected):
    __slots__ = ['annotated_seq']

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return [None] * self.variant_count


@dataclass(frozen=True)
class PamProtCDSTargeton(CDSTargeton[PamVariant, GenomicRange], PamProtected):
    __slots__ = ['annotated_seq', '_codon_to_sgrna_id']

    _SGRNA_IDS: ClassVar[str] = 'sgrna_ids'

    def _get_codon_to_sgrna_id(self) -> Dict[int, str]:
        """Create mapping of codon indices to sgRNA identifiers"""

        return {
            codon_index: variant.sgrna_id
            for variant, codon_index in zip(
                self.annotated_seq.variants,
                self.annotated_seq.get_variant_codon_indices())
        } if self.has_pam_variants else {}

    def __post_init__(self) -> None:
        if self.annotated_seq.contains_same_codon_variants:
            self.annotated_seq.log_same_codon_variants()
            raise ValueError("Multiple PAM protection variants in a single codon!")

        object.__setattr__(self, '_codon_to_sgrna_id', self._get_codon_to_sgrna_id())

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return self.annotated_seq.get_variant_mutation_types(
            codon_table, no_duplicate_codons=True)

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return frozenset({
            self._codon_to_sgrna_id[codon_index]
            for codon_index in self.get_codon_indices(spr)
            if codon_index in self._codon_to_sgrna_id
        })

    @property
    def has_pam_variants(self) -> bool:
        return self.variant_count > 0

    @classmethod
    def from_annotated_sequence(cls, annotated_sequence: CDSAnnotatedSequencePair) -> PamProtCDSTargeton:
        return cls(annotated_sequence, codon_to_sgrna_id)

    def set_metadata_sgrna_ids(self, mutations: MutationCollection) -> None:
        return set_metadata_sgrna_ids(self.frame, self._codon_to_sgrna_id, mutations.df)

    def _assign_pam_sgrna_ids_to_mutations(self, sgrna_ids: FrozenSet[str], mutations: MutationCollection) -> None:
        if not self.has_pam_variants or mutations.is_empty:
            return

        self.set_metadata_sgrna_ids(mutations)

    def _get_sgrna_ids(self, d: Dict[str, Any]) -> FrozenSet[str]:
        sgrna_ids: Optional[FrozenSet[str]] = d.get(self._SGRNA_IDS, None)
        if sgrna_ids is None:
            raise RuntimeError("sgRNA ID's required!")
        if not isinstance(sgrna_ids, frozenset):
            raise TypeError("sgRNA ID's: not a FrozenSet!")
        return sgrna_ids

    def compute_mutations(
        self,
        mutators: FrozenSet[TargetonMutator],
        aux: AuxiliaryTables,
        **kwargs
    ) -> Dict[TargetonMutator, MutationCollection]:
        f"""
        Compute all mutations

        sgRNA ID's should be provided unless no PAM protection was applied.

        Keyword arguments:
        - {self._SGRNA_IDS}: FrozenSet[str]
        """

        if self.has_pam_variants:
            sgrna_ids = self._get_sgrna_ids(kwargs)

        mutations = super().compute_mutations(mutators, aux)

        if self.has_pam_variants:
            f = partial(self._assign_pam_sgrna_ids_to_mutations, sgrna_ids)
            for mc in mutations.values():
                f(mc)

        return mutations
