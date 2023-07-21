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
from typing import Dict, Iterable, List, Optional, Set, FrozenSet, Tuple
import numpy as np
import pandas as pd
from ..metadata_utils import get_mave_nt_pam_from_row, get_mave_nt_ref_from_row, set_pam_extended_ref_alt, set_ref, set_ref_meta
from .codon_table import CodonTable
from .oligo_segment import OligoSegment, TargetonOligoSegment
from .refseq_ranges import ReferenceSequenceRanges
from .targeton import PamProtCDSTargeton
from .base import GenomicRange, TranscriptInfo
from .custom_variants import CustomVariantMutation, CustomVariantMutationCollection, CustomVariantOligoRenderer
from .mutated_sequences import MutationCollection
from .oligo_renderer import BaseOligoRenderer
from .options import Options
from .pam_protected_reference_sequence import PamProtectedReferenceSequence
from .snv_table import AuxiliaryTables
from .variant import CustomVariant
from ..constants import (
    ARRAY_SEPARATOR,
    CUSTOM_MUTATOR,
    META_MAVE_NT,
    META_MAVE_NT_REF,
    META_OLIGO_NAME,
    META_PAM_MUT_ANNOT,
    META_VAR_TYPE,
    META_MUT_POSITION,
    META_REF,
    META_NEW,
    META_MUTATOR,
    META_MSEQ,
    META_OLIGO_LENGTH,
    META_VCF_ALIAS,
    MUTATION_TYPE_NON_CDS
)
from ..enums import MutationType, TargetonMutator
from ..utils import get_constant_category, group_consecutive

MUTATION_TYPE_LABELS: Dict[int, str] = {
    MutationType.SYNONYMOUS.value: 'syn',
    MutationType.MISSENSE.value: 'mis',
    MutationType.NONSENSE.value: 'non'
}

MUTATION_TYPE_CATEGORIES = sorted(MUTATION_TYPE_LABELS.values())
MUTATION_TYPE_CATEGORIES_T = tuple(MUTATION_TYPE_CATEGORIES)


def _decode_mut_type(x) -> str:
    return MUTATION_TYPE_LABELS[x] if not pd.isnull(x) else x


def decode_mut_types(mut_type: pd.Series) -> pd.Series:
    return mut_type.apply(_decode_mut_type)


def decode_mut_types_cat(mut_type: pd.Series) -> pd.Categorical:
    return pd.Categorical(
        decode_mut_types(mut_type),
        categories=MUTATION_TYPE_CATEGORIES)


def encode_pam_mutation_types(mutation_types: Iterable[Optional[MutationType]]) -> str:
    return ARRAY_SEPARATOR.join([
        MUTATION_TYPE_LABELS[x.value] if x is not None else MUTATION_TYPE_NON_CDS
        for x in mutation_types
    ])


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

    def _get_mutated_sequence(self, mseq: str) -> str:
        return f"{self.prefix}{mseq}{self.suffix}"


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

        df: pd.DataFrame = self.mutation_collection.df
        df[META_MUTATOR] = get_constant_category(self.mutator.value, df.shape[0])
        df[META_MUT_POSITION] += self.target_region_start
        df[META_OLIGO_NAME] = self.renderer.get_oligo_names_from_dataframe(df)

        return self.renderer.get_metadata_table(df, options)


@dataclass
class OligoTemplate:
    __slots__ = {
        'ref_ranges',
        'transcript_info',
        'ref_seq',
        'sgrna_ids',
        'custom_variants',
        'adaptor_5',
        'adaptor_3',
        'segments'
    }

    ref_ranges: ReferenceSequenceRanges
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
    def ref_start(self) -> int:
        return self.ref_range.start

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
        # TODO: check the fact these could be the constant regions has no consequence;
        # probably not because they have no mutators assigned.
        return [
            (i, s)
            for i, s in enumerate(self.segments)
            if isinstance(s, TargetonOligoSegment)
        ]

    @property
    def ref_segments(self) -> List[OligoSegment]:
        return self.segments

    def _get_concat_cds_targetons(self) -> List[PamProtCDSTargeton]:
        # Assumption: segments are consecutive (concatenation would fail otherwise)
        targetons = [
            (i, segment.targeton)
            for i, segment in self.target_segments
            if isinstance(segment.targeton, PamProtCDSTargeton)
        ]

        n: int = len(targetons)

        if n == 0:
            return []
        if n == 1:
            return [targetons[0][1]]

        grouped_targetons = group_consecutive(targetons)

        if (
            logging.root.level <= logging.DEBUG
            and any(len(g) > 1 for g in grouped_targetons)
        ):
            logging.debug("Consecutive CDS targeton concatenation required.")

        # Targetons will be extracted as they are from single-targeton groups
        return list(map(PamProtCDSTargeton.concat, grouped_targetons))

    def _get_custom_variant_sgrna_ids(self, variant: CustomVariant) -> FrozenSet[str]:
        spr = variant.get_ref_range(self.strand)
        return frozenset().union(*[
            segment.get_sgrna_ids(spr)
            for segment in self.segments
        ])

    def _get_custom_variant_mutation(self, variant: CustomVariant) -> CustomVariantMutation:
        return CustomVariantMutation(
            variant,
            self.ref_seq.get_variant_corrected_ref(
                variant.base_variant),
            self.ref_seq.apply_variant(
                variant.base_variant, ref_check=False),
            self.ref_ranges.is_range_in_constant_region(
                variant.get_ref_range(self.strand)),
            self._get_custom_variant_sgrna_ids(variant))

    def _compute_custom_variants(self) -> CustomVariantMutationCollection:
        return CustomVariantMutationCollection.from_variants(
            list(map(self._get_custom_variant_mutation, self.custom_variants)))

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
                r[META_VCF_ALIAS],
                r[META_VAR_TYPE],
                r[META_MUT_POSITION],
                r[META_REF] if not pd.isnull(r[META_REF]) else None,
                r[META_NEW] if not pd.isnull(r[META_NEW]) else None), axis=1), dtype='string')
        df['mutator'] = get_constant_category(CUSTOM_MUTATOR, df.shape[0])
        return renderer.get_metadata_table(df, options)

    def _get_pam_variant_annotations(self, codon_table: CodonTable, n: int) -> pd.Categorical:
        s: str = encode_pam_mutation_types(chain.from_iterable([
            segment.get_pam_annotations(codon_table)
            for segment in self.segments
        ]))
        return get_constant_category(s, n, categories=[s])

    def get_mutation_table(self, aux: AuxiliaryTables, options: Options) -> Optional[pd.DataFrame]:
        all_mutations: Optional[pd.DataFrame]

        # Compute mutations per region
        region_mutations: Optional[pd.DataFrame] = pd.concat([
            pd.concat([
                self._get_mutation_collection(
                    i, segment, mutator, mutation_collection).get_metadata_table(options)
                for mutator, mutation_collection in segment.compute_mutations(
                    aux, sgrna_ids=self.sgrna_ids).items()
            ])
            for i, segment in self.target_segments
        ], ignore_index=True) if self.target_segments else None

        # Compute global mutations (custom variants)
        if self.custom_variants:
            global_mutations: pd.DataFrame = self._get_custom_variant_collection(options)
            if global_mutations.shape[0] != len(self.custom_variants):
                raise RuntimeError("Unexpected number of custom variants!")

            # Concatenate to region mutations, if any
            all_mutations = (
                pd.concat([
                    region_mutations,
                    global_mutations
                ]) if region_mutations is not None else
                global_mutations
            )

        else:
            # BEWARE: if empty, the table won't have the expected dtypes
            all_mutations = region_mutations

        if all_mutations is None or all_mutations.empty:
            return None

        # Decode mutation type
        if 'mut_type' in all_mutations.columns:
            # TODO: improve filter performance
            all_mutations.mut_type = decode_mut_types(all_mutations.mut_type)

        # Compute oligonucleotide lengths
        all_mutations[META_OLIGO_LENGTH] = all_mutations[META_MSEQ].str.len().astype(np.int32)

        # Add PAM variant annotations
        rown: int = all_mutations.shape[0]
        all_mutations[META_PAM_MUT_ANNOT] = self._get_pam_variant_annotations(
            aux.codon_table, rown)

        assert 'pam_seq' in all_mutations

        # Remove duplications from the index
        all_mutations = all_mutations.reset_index(drop=True)

        # Concatenate consecutive CDS targetons and discard those unaffected by PAM variants
        pam_prot_cds_targetons = [
            targeton
            for targeton in self._get_concat_cds_targetons()
            if targeton.variant_count > 0
        ]

        set_ref_meta(all_mutations)
        set_ref(all_mutations)
        set_pam_extended_ref_alt(all_mutations, pam_prot_cds_targetons)

        # Add mutation MAVE-HGVS code (REF without PAM protection)
        all_mutations[META_MAVE_NT_REF] = pd.Series(all_mutations.apply(get_mave_nt_ref_from_row, axis=1), dtype='string')

        # Add mutation MAVE-HGVS code (REF and ALT extended to include the PAM codon)
        all_mutations[META_MAVE_NT] = pd.Series(all_mutations.apply(get_mave_nt_pam_from_row, axis=1), dtype='string')

        return all_mutations
