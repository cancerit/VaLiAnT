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
from collections.abc import Sized
from dataclasses import dataclass
from typing import Callable, ClassVar, Dict, List, FrozenSet, Optional
import numpy as np
import pandas as pd
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
from .base import StrandedPositionRange, GenomicRange
from .cdna import CDNA
from .pam_protection import PamProtectedReferenceSequence
from .sequences import Sequence
from .snv_table import AuxiliaryTables
from ..enums import TargetonMutator, VariantType
from ..string_mutators import delete_non_overlapping_3_offset, replace_codons_const
from ..utils import get_constant_category, get_out_of_frame_offset


def get_snv_mutations(sequence: str) -> MutationCollection:
    return MutationCollection.from_variants(
        SingleNucleotideMutatedSequence.from_sequence(
            sequence))


class BaseTargeton(abc.ABC, Sized):
    MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = GENERIC_MUTATORS

    def __init__(self, seq: str, pam_seq: str, pr: StrandedPositionRange) -> None:
        self.seq = seq
        self.pam_seq = pam_seq
        self.pos_range = pr

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def sequence(self) -> str:
        return self.pam_seq

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


class Targeton(BaseTargeton):
    __slots__ = {'ref_sequence'}

    @classmethod
    def from_pam_seq(cls, ref_sequence: PamProtectedReferenceSequence) -> Targeton:
        return cls(
            ref_sequence.sequence,
            ref_sequence.pam_protected_sequence,
            ref_sequence.genomic_range)

    def compute_mutations(self, mutators: FrozenSet[TargetonMutator]) -> Dict[TargetonMutator, MutationCollection]:
        return super()._compute_mutations(mutators)


class CDSTargeton(BaseTargeton):
    __slots__ = {'ref_sequence', 'cds_prefix', 'cds_suffix'}

    MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = GENERIC_MUTATORS | CDS_ONLY_MUTATORS

    SNVRE_MUTATORS: ClassVar[FrozenSet[TargetonMutator]] = frozenset([
        TargetonMutator.SNV_RE
    ])

    def __init__(self, seq: str, pam_seq: str, pr: StrandedPositionRange, cds_prefix: str, cds_suffix: str) -> None:
        super().__init__(seq, pam_seq, pr)
        self.cds_prefix: str = cds_prefix
        self.cds_suffix: str = cds_suffix

    def __post_init__(self) -> None:
        if (len(self) + len(self.cds_prefix) + len(self.cds_suffix)) % 3 != 0:
            raise ValueError("Invalid length for in-frame sequence!")

    @classmethod
    def from_pam_seq(cls, ref_sequence: PamProtectedReferenceSequence, cds_prefix: str, cds_suffix: str) -> CDSTargeton:
        return cls(
            ref_sequence.sequence,
            ref_sequence.pam_protected_sequence,
            ref_sequence.genomic_range,
            cds_prefix,
            cds_suffix)

    @property
    def strand(self) -> str:
        return self.pos_range.strand

    @property
    def cds_sequence(self) -> str:
        return f"{self.cds_prefix}{self.sequence}{self.cds_suffix}"

    @property
    def frame(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

    @property
    def start(self) -> int:
        return self.pos_range.start

    @property
    def cds_sequence_start(self) -> int:
        return self.start - self.frame

    def _add_snv_metadata(
        self,
        aux: AuxiliaryTables,
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
            snvs: MutationCollection = self._add_snv_metadata(aux, snv_meta_full, mutations[TargetonMutator.SNV])
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
