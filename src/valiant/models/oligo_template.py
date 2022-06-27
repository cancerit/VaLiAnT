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
from functools import partial
from typing import Callable, Dict, List, Optional, Set, FrozenSet, Tuple
import numpy as np
import pandas as pd
from .base import GenomicRange, TranscriptInfo
from .custom_variants import CustomVariantMutation, CustomVariantMutationCollection, CustomVariantOligoRenderer
from .mutated_sequences import MutatedSequence, MutationCollection
from .oligo_renderer import BaseOligoRenderer
from .options import Options
from .pam_protection import PamProtectedReferenceSequence
from .snv_table import AuxiliaryTables
from .targeton import BaseTargeton, CDSTargeton, Targeton
from .variant import CustomVariant
from ..constants import (
    CUSTOM_MUTATOR,
    META_OLIGO_NAME,
    META_VAR_TYPE,
    META_MUT_POSITION,
    META_REF,
    META_NEW,
    META_MUTATOR,
    META_MSEQ,
    META_OLIGO_LENGTH
)
from ..enums import MutationType, TargetonMutator
from ..utils import get_constant_category

MUTATION_TYPE_LABELS: Dict[int, str] = {
    MutationType.SYNONYMOUS.value: 'syn',
    MutationType.MISSENSE.value: 'mis',
    MutationType.NONSENSE.value: 'non'
}

MUTATION_TYPE_CATEGORIES = sorted(MUTATION_TYPE_LABELS.values())
MUTATION_TYPE_CATEGORIES_T = tuple(MUTATION_TYPE_CATEGORIES)

EMPTY_MUTATION_TABLE_FIELDS = [
    META_OLIGO_NAME,
    META_VAR_TYPE,
    META_MUT_POSITION,
    META_REF,
    META_NEW,
    META_MUTATOR,
    META_MSEQ,
    META_OLIGO_LENGTH
]


def _decode_mut_type(x) -> str:
    return MUTATION_TYPE_LABELS[x] if not pd.isnull(x) else x


def decode_mut_types(mut_type: pd.Series) -> pd.Series:
    return mut_type.apply(_decode_mut_type)


def decode_mut_types_cat(mut_type: pd.Series) -> pd.Categorical:
    return pd.Categorical(
        decode_mut_types(mut_type),
        categories=MUTATION_TYPE_CATEGORIES)


def get_empty_mutation_table() -> pd.DataFrame:
    return pd.DataFrame(columns=EMPTY_MUTATION_TABLE_FIELDS)


class OligoSegment(abc.ABC):

    @property
    @abc.abstractmethod
    def sequence(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def genomic_range(self) -> GenomicRange:
        pass

    @property
    @abc.abstractmethod
    def pam_protected_sequence(self) -> str:
        pass

    @property
    def start(self) -> int:
        return self.genomic_range.start


@dataclass
class InvariantOligoSegment(OligoSegment):
    __slots__ = {'ref_sequence'}

    ref_sequence: PamProtectedReferenceSequence

    @property
    def sequence(self) -> str:
        return self.ref_sequence.sequence

    @property
    def genomic_range(self) -> GenomicRange:
        return self.ref_sequence.genomic_range

    @property
    def pam_protected_sequence(self) -> str:
        return self.ref_sequence.pam_protected_sequence

    @property
    def ref_seq(self) -> PamProtectedReferenceSequence:
        return self.ref_sequence


@dataclass
class TargetonOligoSegment(OligoSegment):
    __slots__ = {'targeton', 'mutator'}

    targeton: BaseTargeton
    mutators: FrozenSet[TargetonMutator]

    @property
    def sequence(self) -> str:
        return self.targeton.seq

    @property
    def pam_protected_sequence(self) -> str:
        return self.targeton.pam_seq

    @property
    def genomic_range(self) -> GenomicRange:
        # This property can only be used for SGE libraries
        return self.targeton.pos_range  # type: ignore

    def compute_mutations(self, aux: AuxiliaryTables) -> Dict[TargetonMutator, MutationCollection]:
        if isinstance(self.targeton, CDSTargeton):
            return self.targeton.compute_mutations(self.mutators, aux)
        if isinstance(self.targeton, Targeton):
            return self.targeton.compute_mutations(self.mutators)
        raise TypeError("Invalid targeton type!")

    @property
    def start(self) -> int:
        return self.targeton.start


@dataclass(init=False)
class RegionOligoRenderer(BaseOligoRenderer):
    __slots__ = {
        'ref_seq',
        'gene_id',
        'transcript_id',
        'adaptor_5',
        'adaptor_3',
        'prefix',
        'suffix',
        '_oligo_name_prefix'
    }

    prefix: str
    suffix: str

    def __init__(
        self,
        ref_seq: PamProtectedReferenceSequence,
        gene_id: str,
        transcript_id: str,
        adaptor_5: str,
        adaptor_3: str,
        prefix: str,
        suffix: str
    ) -> None:
        super().__init__(ref_seq, gene_id, transcript_id, adaptor_5, adaptor_3)
        self.prefix = prefix
        self.suffix = suffix

    def get_oligo_name(self, mutator: TargetonMutator, start_offset: int, ms: MutatedSequence) -> str:
        return super()._get_oligo_name(
            ms.type.value, mutator.value, start_offset + ms.position, ms.ref, ms.new)

    def _render_mutated_ref_sequence(self, mseq: str) -> str:
        return f"{self.prefix}{mseq}{self.suffix}"

    def _render_mutated_sequence(self, mseq: str) -> str:
        return super()._render_mutated_sequence(
            self._render_mutated_ref_sequence(mseq))

    def _render_mutated_sequence_rc(self, mseq: str) -> str:
        return super()._render_mutated_sequence_rc(
            self._render_mutated_ref_sequence(mseq))


@dataclass
class OligoMutationCollection:
    __slots__ = {
        'renderer',
        'mutator',
        'mutation_collection',
        'target_region_start'
    }

    renderer: RegionOligoRenderer
    mutator: TargetonMutator
    mutation_collection: MutationCollection
    target_region_start: int

    def get_metadata_table(self, options: Options) -> pd.DataFrame:
        if self.mutation_collection.is_empty:
            raise RuntimeError(
                f"Empty mutation collection for mutator '{self.mutator}'!")

        get_oligo_name: Callable[[MutatedSequence], str] = partial(
            self.renderer.get_oligo_name,
            self.mutator,
            self.target_region_start)

        df: pd.DataFrame = self.mutation_collection.df
        df['oligo_name'] = pd.Series(
            map(get_oligo_name, self.mutation_collection.mutations),
            dtype='string')
        df['mutator'] = get_constant_category(self.mutator.value, df.shape[0])
        df.mut_position += self.target_region_start

        return self.renderer.get_metadata_table(df, options)


@dataclass
class OligoTemplate:
    __slots__ = {
        'transcript_info',
        'ref_seq',
        'sgrna_ids',
        'custom_variants',
        'adaptor_5',
        'adaptor_3',
        'segments'
    }

    transcript_info: Optional[TranscriptInfo]
    ref_seq: PamProtectedReferenceSequence
    sgrna_ids: FrozenSet[str]
    custom_variants: Set[CustomVariant]
    adaptor_5: Optional[str]
    adaptor_3: Optional[str]
    segments: List[OligoSegment]

    @property
    def transcript_id(self) -> Optional[str]:
        return self.transcript_info.transcript_id if self.transcript_info else None

    @property
    def gene_id(self) -> Optional[str]:
        return self.transcript_info.gene_id if self.transcript_info else None

    @property
    def ref_range(self) -> GenomicRange:
        return self.ref_seq.genomic_range

    @property
    def strand(self) -> str:
        return self.ref_range.strand

    @property
    def name(self) -> str:
        return '_'.join([
            self.ref_range.chromosome,
            str(self.ref_range.start),
            str(self.ref_range.end),
            'plus' if self.ref_range.strand == '+' else 'minus'
        ])

    @property
    def segment_strings(self) -> List[str]:
        """Render PAM protected sequences to populate the oligonucleotide prefix and suffix"""

        return [s.pam_protected_sequence for s in self.segments]

    @property
    def target_segments(self) -> List[Tuple[int, TargetonOligoSegment]]:
        return [
            (i, s)
            for i, s in enumerate(self.segments)
            if isinstance(s, TargetonOligoSegment)
        ]

    @property
    def ref_segments(self) -> List[OligoSegment]:
        return self.segments

    def _compute_custom_variants(self) -> CustomVariantMutationCollection:
        return CustomVariantMutationCollection.from_variants([
            CustomVariantMutation(
                variant, self.ref_seq.apply_variant(
                    variant.base_variant, ref_check=False))
            for variant in self.custom_variants
        ])

    def _get_mutation_collection(
        self,
        i: int,
        target_segment: TargetonOligoSegment,
        mutator: TargetonMutator,
        mutation_collection: MutationCollection
    ) -> OligoMutationCollection:
        segment_strings: List[str] = self.segment_strings
        prefix: str = ''.join(segment_strings[:i])
        suffix: str = ''.join(segment_strings[i + 1:])

        renderer: RegionOligoRenderer = RegionOligoRenderer(
            self.ref_seq,
            self.gene_id or '',
            self.transcript_id or '',
            self.adaptor_5 or '',
            self.adaptor_3 or '',
            prefix,
            suffix)

        return OligoMutationCollection(
            renderer,
            mutator,
            mutation_collection,
            target_segment.start)

    def _get_custom_variant_collection(self, options: Options) -> pd.DataFrame:
        renderer: CustomVariantOligoRenderer = CustomVariantOligoRenderer(
            self.ref_seq,
            self.gene_id or '',
            self.transcript_id or '',
            self.adaptor_5 or '',
            self.adaptor_3 or '')

        mc: CustomVariantMutationCollection = self._compute_custom_variants()

        if mc.is_empty:
            return pd.DataFrame()

        df: pd.DataFrame = mc.df
        df['oligo_name'] = pd.Series(
            df.apply(lambda r: renderer.get_oligo_name(
                r.vcf_alias,
                r.var_type,
                r.mut_position,
                r.ref if not pd.isnull(r.ref) else None,
                r.new if not pd.isnull(r.new) else None), axis=1), dtype='string')
        df['mutator'] = get_constant_category(CUSTOM_MUTATOR, df.shape[0])
        return renderer.get_metadata_table(df, options)

    def get_mutation_table(self, aux: AuxiliaryTables, options: Options) -> pd.DataFrame:
        if not self.target_segments:
            return get_empty_mutation_table()

        # Compute mutations per region
        region_mutations: pd.DataFrame = pd.concat([
            pd.concat([
                self._get_mutation_collection(
                    i, segment, mutator, mutation_collection).get_metadata_table(options)
                for mutator, mutation_collection in segment.compute_mutations(aux).items()
            ])
            for i, segment in self.target_segments
        ], ignore_index=True)

        # Compute global mutations (custom variants)
        if self.custom_variants:
            global_mutations: pd.DataFrame = self._get_custom_variant_collection(options)
            if global_mutations.shape[0] != len(self.custom_variants):
                raise RuntimeError("Unexpected number of custom variants!")
            all_mutations = pd.concat([region_mutations, global_mutations])
        else:
            all_mutations = region_mutations

        # Decode mutation type
        if 'mut_type' in all_mutations.columns:
            # TODO: improve filter performance
            all_mutations.mut_type = decode_mut_types(all_mutations.mut_type)

        # Compute oligonucleotide lengths
        all_mutations['oligo_length'] = all_mutations.mseq.str.len().astype(np.int32)

        return all_mutations
