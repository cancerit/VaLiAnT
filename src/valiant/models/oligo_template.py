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
from functools import partial
from typing import Callable, Dict, List, Optional, Set, Tuple
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
from ..constants import CUSTOM_MUTATOR
from ..enums import MutationType, TargetonMutator
from ..utils import get_constant_category

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
    mutators: Set[TargetonMutator]

    @property
    def sequence(self) -> str:
        return self.targeton.seq

    @property
    def pam_protected_sequence(self) -> str:
        return self.targeton.pam_seq

    @property
    def genomic_range(self) -> GenomicRange:
        return self.targeton.pos_range

    def compute_mutations(self, aux: AuxiliaryTables) -> Dict[TargetonMutator, MutationCollection]:
        if isinstance(self.targeton, CDSTargeton):
            return self.targeton.compute_mutations(self.mutators, aux)
        if isinstance(self.targeton, Targeton):
            return self.targeton.compute_mutations(self.mutators)
        raise TypeError("Invalid targeton type!")

    @property
    def start(self) -> int:
        return self.targeton.pos_range.start


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
    sgrna_ids: Set[str]
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
            return pd.DataFrame(columns=[
                'oligo_name',
                'var_type',
                'mut_position',
                'ref',
                'new',
                'mutator',
                'mseq',
                'oligo_length'
            ])

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
