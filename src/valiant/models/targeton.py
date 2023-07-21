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
from typing import Any, Callable, ClassVar, Dict, Generic, List, FrozenSet, Optional, Tuple, Type, TypeVar
import numpy as np
import pandas as pd

from ..constants import GENERIC_MUTATORS, CDS_ONLY_MUTATORS
from ..enums import MutationType, TargetonMutator, VariantType
from ..sgrna_utils import set_metadata_sgrna_ids, set_metadata_sgrna_ids_empty
from ..string_mutators import delete_non_overlapping_3_offset, replace_codons_const
from ..utils import get_constant_category, get_out_of_frame_offset
from .annotated_sequence import AnnotatedSequencePair, AnnotatedSequenceT, CDSAnnotatedSequencePair, RangeT, VariantT
from .base import GenomicRange, StrandedPositionRange
from .codon_table import CodonTable, STOP_CODE
from .pam_protected_reference_sequence import PamProtectedReferenceSequence
from .pam_protection import PamVariant, get_position_to_sgrna_ids
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
    def variants(self) -> List[VariantT]:
        return self.annotated_seq.variants

    @property
    def variant_count(self) -> int:
        return self.annotated_seq.variant_count

    def _get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return frozenset({
            variant.sgrna_id
            for variant in self.annotated_seq.get_variants_in_range(spr)
        })

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
        return super()._get_sgrna_ids(spr)

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
        return self.annotated_seq.start

    @property
    def end(self) -> int:
        return self.annotated_seq.end

    @property
    def frame(self) -> int:
        return self.annotated_seq.frame

    def get_codon_indices(self, spr: StrandedPositionRange) -> List[int]:
        return self.annotated_seq.get_codon_indices_in_range(spr)

    def get_liminal_codon_indices(self, spr: StrandedPositionRange) -> Tuple[int, int]:
        codon_indices = self.get_codon_indices(spr)
        return codon_indices[0], codon_indices[-1]

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return super()._get_sgrna_ids(spr)

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

    def _get_rel_position_to_sgrna_ids(self, start: int, variants: List[PamVariant]) -> Dict[int, str]:
        return {
            (k - start): v
            for k, v in get_position_to_sgrna_ids(variants).items()
        }

    @abc.abstractmethod
    def get_position_to_sgrna_ids(self) -> Dict[int, str]:
        pass

    @abc.abstractmethod
    def _set_metadata_sgrna_ids(self, mutations: MutationCollection) -> None:
        pass


@dataclass(frozen=True)
class PamProtTargeton(Targeton[PamVariant, GenomicRange], PamProtected):
    __slots__ = ['annotated_seq']

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return [None] * self.variant_count

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        return super()._get_sgrna_ids(spr)

    def get_position_to_sgrna_ids(self) -> Dict[int, str]:
        return get_position_to_sgrna_ids(self.variants)

    def _set_metadata_sgrna_ids(self, mutations: MutationCollection) -> None:
        if not mutations.is_empty:
            position_to_sgrna_ids = self._get_rel_position_to_sgrna_ids(self.start, self.variants)
            set_metadata_sgrna_ids(0, position_to_sgrna_ids, {}, mutations.df)
        else:
            set_metadata_sgrna_ids_empty(mutations.df)

    def compute_mutations(
        self,
        mutators: FrozenSet[TargetonMutator],
        **kwargs
    ) -> Dict[TargetonMutator, MutationCollection]:
        f"""
        Compute all mutations

        sgRNA ID's should be provided unless no PAM protection was applied.
        """

        mutations = super().compute_mutations(mutators)

        for mc in mutations.values():
            self._set_metadata_sgrna_ids(mc)

        return mutations


@dataclass(frozen=True)
class PamProtCDSTargeton(CDSTargeton[PamVariant, GenomicRange], PamProtected):
    __slots__ = ['annotated_seq', '_codon_to_pam_variant', '_codon_to_pam_codon']

    def get_position_to_sgrna_ids(self) -> Dict[int, str]:
        return get_position_to_sgrna_ids(self.variants)

    @classmethod
    def concat(cls, targetons: List[PamProtCDSTargeton]) -> PamProtCDSTargeton:
        n: int = len(targetons)
        if n == 0:
            raise ValueError("No targetons to merge!")
        if n == 1:
            return targetons[0]

        # Set variables based on the first targeton
        targeton: PamProtCDSTargeton = targetons[0]
        chromosome: str = targeton.pos_range.chromosome
        strand: str = targeton.pos_range.strand
        start: int = targeton.pos_range.start
        cds_prefix: str = targeton.annotated_seq.cds_prefix
        ref_seq: str = targeton.annotated_seq.ref_seq
        alt_seq: str = targeton.annotated_seq.alt_seq
        pam_variants: List[PamVariant] = targeton.annotated_seq.variants

        prev_end: int = targeton.pos_range.end
        prev_suffix: str = targeton.annotated_seq.cds_suffix

        # Set variables based on the last targeton
        targeton = targetons[-1]
        end: int = targeton.pos_range.end
        cds_suffix: str = targeton.annotated_seq.cds_suffix

        pos_range: GenomicRange
        for i in range(1, n):
            targeton = targetons[i]
            pos_range = targeton.pos_range
            if (
                pos_range.chromosome != chromosome
                or pos_range.strand != strand
                or pos_range.start != prev_end + 1
                or targeton.annotated_seq.ref_seq[:len(prev_suffix)] != prev_suffix
            ):
                raise ValueError("Non-consecutive targetons!")

            ref_seq += targeton.annotated_seq.ref_seq
            alt_seq += targeton.annotated_seq.alt_seq
            pam_variants += targeton.annotated_seq.variants

            prev_end = pos_range.end
            prev_suffix = targeton.annotated_seq.cds_suffix

        return PamProtCDSTargeton(CDSAnnotatedSequencePair(
            GenomicRange(chromosome, start, end, strand),
            ref_seq,
            alt_seq,
            pam_variants,
            cds_prefix,
            cds_suffix))

    def _get_codon_to_pam_variant(self, codon_indices: List[int]) -> Dict[int, PamVariant]:
        """Create mapping of codon indices to sgRNA identifiers"""

        return {
            codon_index: variant
            for variant, codon_index in zip(
                self.variants,
                codon_indices)
        }

    def _get_codon_to_pam_codon(self, codon_indices: List[int]) -> Dict[int, str]:
        """Create mapping of codon indices to PAM-protected codon sequences"""

        d = self.annotated_seq.get_indexed_alt_codons()
        return {
            codon_index: d[codon_index]
            for codon_index in codon_indices
        }

    def _setattr(self, attr: str, value: Any) -> None:
        object.__setattr__(self, attr, value)

    def __post_init__(self) -> None:
        if self.annotated_seq.contains_same_codon_variants:
            self.annotated_seq.log_same_codon_variants()
            raise ValueError("Multiple PAM protection variants in a single codon!")

        codon_to_pam_variant: Dict[int, PamVariant] = {}
        codon_to_pam_codon: Dict[int, str] = {}

        if self.has_pam_variants:
            codon_indices = self.annotated_seq.get_variant_codon_indices()
            codon_to_pam_variant = self._get_codon_to_pam_variant(codon_indices)
            codon_to_pam_codon = self._get_codon_to_pam_codon(codon_indices)

        self._setattr('_codon_to_pam_variant', codon_to_pam_variant)
        self._setattr('_codon_to_pam_codon', codon_to_pam_codon)

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return self.annotated_seq.get_variant_mutation_types(
            codon_table, no_duplicate_codons=True)

    def _get_pam_codon_length(self, codon_index: int) -> int:
        return len(self._codon_to_pam_codon[codon_index])

    def _get_pam_codon_theoretical_range(self, codon_index: int) -> Tuple[int, int]:
        start: int = (codon_index * 3) + self.pos_range.start - self.frame
        end: int = start + 2
        return start, end

    def _get_pam_codon_position(self, codon_index: int, offset_f: Callable[[int], int]) -> int:
        return (codon_index * 3) + self.pos_range.start - self.frame + offset_f(
            self._get_pam_codon_length(codon_index))

    def _get_pam_codon_start(self, codon_index: int) -> int:
        return self._get_pam_codon_position(
            codon_index, lambda codon_length: 3 - codon_length)

    def _get_pam_codon_end(self, codon_index: int) -> int:
        return self._get_pam_codon_position(
            codon_index, lambda codon_length: codon_length - 1)

    def _is_codon_pam_protected(self, codon_index: int) -> bool:
        return codon_index in self._codon_to_pam_variant

    def _get_pam_variant_position_by_codon_index(self, codon_index: int) -> int:
        return self._codon_to_pam_variant[codon_index].genomic_position.position

    def get_pam_ext_start(self, pos: int) -> int:
        """Get genomic position at the start of the PAM-protected codon"""

        codon_index: int = self.annotated_seq.get_codon_index(pos)

        return (
            min(
                pos,
                self._get_pam_variant_position_by_codon_index(codon_index)
            ) if self._is_codon_pam_protected(codon_index) else
            pos
        )

    def get_pam_ext_ref_end(self, pos: int) -> int:
        """Get genomic position at the end of the PAM-protected codon"""

        codon_index: int = self.annotated_seq.get_codon_index(pos)

        return (
            max(
                pos,
                self._get_pam_variant_position_by_codon_index(codon_index)
            ) if self._is_codon_pam_protected(codon_index) else
            pos
        )

    def get_sgrna_ids(self, spr: StrandedPositionRange) -> FrozenSet[str]:
        cds = frozenset({
            self._codon_to_pam_variant[codon_index].sgrna_id
            for codon_index in self.get_codon_indices(spr)
            if codon_index in self._codon_to_pam_variant
        }) if self._codon_to_pam_variant else frozenset()
        ncd = super()._get_sgrna_ids(spr)
        return frozenset.union(cds, ncd)

    @property
    def has_pam_variants(self) -> bool:
        return self.variant_count > 0

    def _set_metadata_sgrna_ids(self, mutations: MutationCollection) -> None:
        if not mutations.is_empty:
            codon_to_sgrna_ids: Dict[int, str] = {
                codon_index: variant.sgrna_id
                for codon_index, variant in self._codon_to_pam_variant.items()
            } if self.has_pam_variants else {}
            position_to_sgrna_ids = self._get_rel_position_to_sgrna_ids(self.start, self.variants)
            set_metadata_sgrna_ids(self.frame, position_to_sgrna_ids, codon_to_sgrna_ids, mutations.df)
        else:
            set_metadata_sgrna_ids_empty(mutations.df)

    def compute_mutations(
        self,
        mutators: FrozenSet[TargetonMutator],
        aux: AuxiliaryTables,
        **kwargs
    ) -> Dict[TargetonMutator, MutationCollection]:
        """
        Compute all mutations

        sgRNA ID's should be provided unless no PAM protection was applied.
        """

        mutations = super().compute_mutations(mutators, aux)

        for mc in mutations.values():
            self._set_metadata_sgrna_ids(mc)

        return mutations
