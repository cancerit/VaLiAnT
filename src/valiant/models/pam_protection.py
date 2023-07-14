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
import logging
from typing import Dict, Iterable, List, Optional, Set, FrozenSet
import pandas as pd
from pyranges import PyRanges
from ..loaders.vcf import get_vcf
from ..utils import get_id_column
from .base import GenomicRange
from .sequences import ReferenceSequence
from .variant import BaseVariant, get_variant, SubstitutionVariant, apply_variants


@dataclass(frozen=True)
class PamProtectedReferenceSequence(ReferenceSequence):
    __slots__ = [
        'sequence',
        'genomic_range',
        'pam_protected_sequence',
        'bg_sequence',
        'pam_variants',
        'bg_variants'
    ]

    pam_protected_sequence: str
    bg_sequence: str
    pam_variants: List[PamVariant]
    bg_variants: List[BaseVariant]

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.pam_protected_sequence) != len(self.genomic_range):
            raise ValueError("PAM protected sequence and genomic range have different lengths!")

    @classmethod
    def from_reference_sequence(
        cls,
        ref_seq: ReferenceSequence,
        pam_variants_: Iterable[PamVariant],
        bg_variants_: Iterable[BaseVariant] = list()
    ) -> PamProtectedReferenceSequence:
        pam_variants: List[PamVariant] = sorted(
            pam_variants_, key=lambda x: x.genomic_position.position)
        bg_variants: List[BaseVariant] = sorted(
            bg_variants_, key=lambda x: x.genomic_position.position)

        # Apply background variants
        bg_seq: str = apply_variants(ref_seq, bg_variants, ref_check=True)
        if len(bg_seq) != len(ref_seq):
            raise NotImplementedError("Background variants altering the length of the reference sequence!")

        # Apply PAM sequence variants
        pam_seq: str = apply_variants(
            ReferenceSequence(bg_seq, ref_seq.genomic_range),
            pam_variants,
            ref_check=not bool(bg_variants))

        return cls(
            ref_seq.sequence,
            ref_seq.genomic_range,
            pam_seq,
            bg_seq,
            pam_variants,
            bg_variants)

    def get_subsequence(self, genomic_range: GenomicRange) -> PamProtectedReferenceSequence:
        ref_seq: ReferenceSequence = super().get_subsequence(genomic_range)
        pam_variants = [
            variant
            for variant in self.pam_variants
            if genomic_range.contains_position(variant.genomic_position)
        ]
        bg_variants = [
            variant
            for variant in self.bg_variants
            if genomic_range.contains_position(variant.genomic_position)
        ]
        return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_variants, bg_variants)

    def _validate_variant_position(self, variant: BaseVariant) -> None:
        if not self.genomic_range.contains_position(variant.genomic_position):
            raise ValueError("Variant not in genomic range!")

    def apply_variant(self, variant: BaseVariant, ref_check: bool = False) -> str:
        self._validate_variant_position(variant)
        return variant.mutate(self.pam_protected_sequence, self.genomic_range.start, ref_check=ref_check)

    def get_variant_corrected_ref(self, variant: BaseVariant) -> Optional[str]:
        self._validate_variant_position(variant)
        return variant.get_corrected_ref(self.pam_protected_sequence, self.genomic_range.start)


@dataclass(frozen=True)
class PamVariant(SubstitutionVariant):
    __slots__ = {'genomic_position', 'ref', 'alt', 'sgrna_id'}

    sgrna_id: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.ref) != 1 or len(self.alt) != 1:
            raise ValueError(
                "Only single-nucleotide substitutions are allowed "
                "for the purposes of PAM protection!")

    @classmethod
    def from_substitution(cls, var: SubstitutionVariant, sgrna_id: str) -> PamVariant:
        if not isinstance(var, SubstitutionVariant):
            raise ValueError("PAM protection variants must be substitutions!")
        return PamVariant(var.genomic_position, var.ref, var.alt, sgrna_id)


class PamProtectionVariantRepository:
    __slots__ = {'sgrna_ids', '_variants', '_ranges'}

    def __init__(self, sgrna_ids: FrozenSet[str] = None) -> None:
        self.sgrna_ids = sgrna_ids or frozenset()
        self._variants: Dict[str, Set[PamVariant]] = {
            sgrna_id: set()
            for sgrna_id in sgrna_ids
        } if sgrna_ids else {}
        self._ranges: Optional[PyRanges] = None

    @property
    def count(self) -> int:
        return sum(map(len, self._variants.values()))

    def get_sgrna_variants(self, sgrna_id: str) -> FrozenSet[PamVariant]:
        try:
            return frozenset(self._variants[sgrna_id])

        except KeyError:
            raise RuntimeError(f"sgRNA ID '{sgrna_id} not found!'")

    def get_sgrna_variants_bulk(self, sgrna_ids: FrozenSet[str]) -> FrozenSet[PamVariant]:
        empty: FrozenSet[PamVariant] = frozenset()
        return empty.union(*[
            self.get_sgrna_variants(sgrna_id)
            for sgrna_id in sgrna_ids
        ])

    def load(self, fp: str) -> None:

        # Load variants from VCF
        with get_vcf(fp) as variant_file:
            for record in variant_file.fetch():
                sgrna_id: str = record.info['SGRNA'].strip()
                if sgrna_id in self._variants:
                    self._variants[sgrna_id].add(
                        PamVariant.from_substitution(get_variant(record), sgrna_id))

        # Log loaded variant statistics
        logging.debug("Collected %d PAM protection variants." % self.count)
        for sgrna_id in self.sgrna_ids:
            if not self._variants[sgrna_id]:
                logging.info("No PAM protection variants for sgRNA %s." % sgrna_id)

        # Populate genomic range table
        df: pd.DataFrame = pd.DataFrame.from_records([
            (*variant.get_pyrange_record(), sgrna_id)
            for sgrna_id, variants in self._variants.items()
            for variant in variants
        ], columns=[
            'Chromosome',
            'Start',
            'End',
            'sgrna_id'
        ])
        df.sgrna_id = df.sgrna_id.astype('category')
        df['variant_id'] = get_id_column(df.shape[0])
        self._ranges = PyRanges(df)


def compute_pam_protected_sequence(
    ref_seq: ReferenceSequence,
    pam_variants: Set[PamVariant],
    bg_variants: Set[BaseVariant]
) -> PamProtectedReferenceSequence:
    return PamProtectedReferenceSequence.from_reference_sequence(
        ref_seq, pam_variants, bg_variants)


def get_position_to_sgrna_ids(pam_variants: List[PamVariant]) -> Dict[int, str]:
    return {
        variant.genomic_position.position: variant.sgrna_id
        for variant in pam_variants
    }
