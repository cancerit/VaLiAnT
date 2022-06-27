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
from typing import Callable, ClassVar, Dict, Generic, List, FrozenSet, Optional
import numpy as np
import pandas as pd

from valiant.models.annotated_sequence import AnnotatedSequencePair, BaseAnnotatedSequencePair, CDSAnnotatedSequencePair, VariantT
from valiant.models.base import StrandedPositionRange
from valiant.models.pam_protection import PamProtectedReferenceSequence, PamVariant
from .codon_table import CodonTable, STOP_CODE
from ..constants import GENERIC_MUTATORS, CDS_ONLY_MUTATORS
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
from ..string_mutators import delete_non_overlapping_3_offset, replace_codons_const
from ..utils import get_constant_category, get_out_of_frame_offset


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
class Targeton(AnnotatedSequencePair, BaseTargeton, Generic[VariantT]):
    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants']

    @property
    def sequence(self) -> str:
        return self.alt_seq

    @classmethod
    def build_without_variants(
        cls,
        pos_range: StrandedPositionRange,
        ref_seq: str
    ) -> Targeton:
        return cls(pos_range, ref_seq, ref_seq, [])

    @classmethod
    def from_pam_seq(
        cls,
        pam_seq: PamProtectedReferenceSequence
    ) -> Targeton:
        return cls(
            pam_seq.genomic_range,
            pam_seq.sequence,
            pam_seq.pam_protected_sequence,
            pam_seq.pam_variants)

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator]) -> Dict[TargetonMutator, MutationCollection]:
        return super()._compute_mutations(mutators)


@dataclass(frozen=True)
class CDSTargeton(CDSAnnotatedSequencePair, BaseTargeton, Generic[VariantT]):
    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants', 'cds_prefix', 'cds_suffix']

    MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = GENERIC_MUTATORS | CDS_ONLY_MUTATORS

    SNVRE_MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = frozenset([
        TargetonMutator.SNV_RE
    ])

    @classmethod
    def build_without_variants(
        cls,
        pos_range: StrandedPositionRange,
        ref_seq: str,
        cds_prefix: str,
        cds_suffix: str
    ) -> CDSTargeton:
        return cls(pos_range, ref_seq, ref_seq, [], cds_prefix, cds_suffix)

    @classmethod
    def from_pam_seq(
        cls,
        pam_seq: PamProtectedReferenceSequence,
        cds_prefix: str,
        cds_suffix: str
    ) -> CDSTargeton:
        return cls(
            pam_seq.genomic_range,
            pam_seq.sequence,
            pam_seq.pam_protected_sequence,
            pam_seq.pam_variants,
            cds_prefix,
            cds_suffix)

    @property
    def sequence(self) -> str:
        return self.alt_seq

    @property
    def strand(self) -> str:
        return self.pos_range.strand

    @property
    def cds_sequence(self) -> str:
        return self.ext_alt_seq

    @property
    def start(self) -> int:
        return self.pos_range.start

    @property
    def cds_sequence_start(self) -> int:
        return self.start - self.frame

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
        return MutationCollection(
            df=snv_joint, mutations=snvs.mutations)

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

        # TODO: avoid generating mutation list when possible
        return MutationCollection(df=df, mutations=[
            SingleCodonMutatedSequence(r.mut_position, r.mseq, r.ref, r.new)
            for r in df.itertuples()
        ])

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
        return MutationCollection(df=df, mutations=[
            SingleCodonMutatedSequence(r.mut_position, r.mseq, r.ref, r.new)
            for r in df.itertuples()
        ])

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator], aux: AuxiliaryTables) -> Dict[TargetonMutator, MutationCollection]:

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
class PamProtTargeton(Targeton[PamVariant], PamProtected):
    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants']

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return [None] * self.variant_count

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator]) -> Dict[TargetonMutator, MutationCollection]:
        return super().compute_mutations(mutators)


@dataclass(frozen=True)
class PamProtCDSTargeton(CDSTargeton[PamVariant], PamProtected):
    __slots__ = ['pos_range', 'ref_seq', 'cds_prefix', 'cds_suffix', '_alt_seq', '_variants']

    def __post_init__(self) -> None:
        if self.contains_same_codon_variants:
            self.log_same_codon_variants()
            raise ValueError("Multiple PAM protection variants in a single codon!")

    def get_pam_variant_annotations(self, codon_table: CodonTable) -> List[Optional[MutationType]]:
        return self.get_variant_mutation_types(
            codon_table, no_duplicate_codons=True)

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator], aux: AuxiliaryTables) -> Dict[TargetonMutator, MutationCollection]:
        return super().compute_mutations(mutators, aux)
