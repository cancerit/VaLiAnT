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

from .background_variants import GenomicPositionOffsets
from ..metadata_utils import get_mave_nt_pam_from_row, get_mave_nt_ref_from_row, set_pam_extended_ref_alt, set_ref, set_ref_meta
from .codon_table import CodonTable
from .oligo_segment import OligoSegment, TargetonOligoSegment
from .refseq_ranges import ReferenceSequenceRanges
from .targeton import PamProtCDSTargeton
from .base import GenomicRange, PositionRange, TranscriptInfo
from .custom_variants import CustomVariantMutation, CustomVariantMutationCollection
from .custom_variant_oligo_renderer import CustomVariantOligoRenderer
from .mutated_sequences import MutationCollection
from .new_pam import PamBgAltSeqBuilder, CdsPamBgAltSeqBuilder, PamBgAltSeqBuilderT
from .oligo_renderer import BaseOligoRenderer
from .options import Options
from .snv_table import AuxiliaryTables
from .custom_variants import CustomVariant
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
from ..utils import get_constant_category, get_nullable_field, group_consecutive, map_filter_log

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
        ref_seq: PamBgAltSeqBuilder,
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

        df: pd.DataFrame = self.mutation_collection.df  # type: ignore
        df[META_MUTATOR] = get_constant_category(self.mutator.value, df.shape[0])
        df[META_MUT_POSITION] += self.target_region_start
        df[META_OLIGO_NAME] = self.renderer.get_oligo_names_from_dataframe(df)

        return self.renderer.get_metadata_table(df, options)


@dataclass
class OligoTemplate:
    __slots__ = {
        'ref_ranges',
        'gpo',
        'transcript_info',
        'ref_seq',
        'sgrna_ids',
        'custom_variants',
        'adaptor_5',
        'adaptor_3',
        'segments'
    }

    ref_ranges: ReferenceSequenceRanges  # corrected for background
    gpo: Optional[GenomicPositionOffsets]
    transcript_info: Optional[TranscriptInfo]
    ref_seq: PamBgAltSeqBuilder
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
        return self.ref_seq.pos_range

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
            self.ref_seq.mutate(
                variant.base_variant),
            self.ref_ranges.is_range_in_constant_region(
                variant.get_ref_range(self.strand)),
            self._get_custom_variant_sgrna_ids(variant))

    def _get_custom_variant_range(self, variant: CustomVariant) -> PositionRange:
        ref_gr = variant.base_variant.ref_range
        if not self.gpo:
            return ref_gr

        alt_gr = self.gpo.ref_to_alt_range(ref_gr)
        logging.debug(f"REF {ref_gr} -> ALT {alt_gr}")
        return PositionRange(alt_gr.start, alt_gr.end)

    def _compute_custom_variants(self) -> CustomVariantMutationCollection:

        def non_overlapping(variant: CustomVariant) -> bool:
            # Convert to ALT coordinates before checking
            return not self.ref_seq.overlaps_bg_range(
                self._get_custom_variant_range(variant))

        def log_overlapping(variant: CustomVariant) -> str:
            return f"Custom variant at {variant.base_variant.genomic_position} " \
                "overlaps with background variant: discarded!"

        return CustomVariantMutationCollection.from_variants(
            map_filter_log(
                items=self.custom_variants,
                map_f=self._get_custom_variant_mutation,
                ft_f=non_overlapping,
                log_f=log_overlapping,
                sort_f=hash,
                level=logging.WARNING))

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

    def _correct_mut_pos(self, df: pd.DataFrame) -> None:
        """Convert the ALT mutation positions to REF (if they can differ)"""

        if self.gpo:
            df[META_MUT_POSITION] = df[META_MUT_POSITION].apply(
                self.gpo.alt_to_ref_position).astype(np.int32)

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

        # Correct the mutation position before generating the name
        self._correct_mut_pos(df)

        df[META_OLIGO_NAME] = pd.Series(
            df.apply(lambda r: renderer.get_oligo_name(
                r[META_VCF_ALIAS],
                r[META_VAR_TYPE],
                r[META_MUT_POSITION],
                get_nullable_field(r, META_REF),
                get_nullable_field(r, META_NEW)), axis=1), dtype='string')
        df[META_MUTATOR] = get_constant_category(CUSTOM_MUTATOR, df.shape[0])
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
            df for _, df in sorted(list(chain.from_iterable([
                [
                    (
                        mutator, self._get_mutation_collection(
                            i, segment, mutator, mutation_collection).get_metadata_table(options)
                    )
                    for mutator, mutation_collection in segment.compute_mutations(
                        aux, sgrna_ids=self.sgrna_ids).items()
                ]
                for i, segment in self.target_segments
            ])), key=lambda t: t[0].value)
        ], ignore_index=True) if self.target_segments else None

        if region_mutations is not None:

            # Correct the mutation position before generating the name
            self._correct_mut_pos(region_mutations)

        # Compute global mutations (custom variants)
        if self.custom_variants:
            global_mutations: pd.DataFrame = self._get_custom_variant_collection(options)

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
