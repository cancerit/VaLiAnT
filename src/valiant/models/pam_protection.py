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
from dataclasses import dataclass
from itertools import chain
import logging
from typing import Dict, Iterable, List, Optional, Set, FrozenSet
import pandas as pd
from pyranges import PyRanges
from pysam import VariantRecord
from .base import GenomicRange
from .sequences import ReferenceSequence
from .variant import BaseVariant, get_variant, SubstitutionVariant
from ..loaders.vcf import get_vcf
from ..utils import get_id_column


def apply_pam_variants(ref_seq: ReferenceSequence, pam_variants: List[PamVariant]) -> str:
    pam_seq: str = ref_seq.sequence

    for variant in pam_variants:

        # Validate variant genomic position relative to the sequence's
        offset: int = variant.get_ref_offset(ref_seq)

        # Update PAM-protected sequence
        pam_seq = variant.mutate_from(pam_seq, offset, ref_check=True)

    return pam_seq


@dataclass(frozen=True)
class PamProtectedReferenceSequence(ReferenceSequence):
    __slots__ = {'sequence', 'genomic_range', 'pam_protected_sequence', 'pam_variants'}

    pam_protected_sequence: str
    pam_variants: List[PamVariant]

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.pam_protected_sequence) != len(self.genomic_range):
            raise ValueError("PAM protected sequence and genomic range have different lengths!")

    @classmethod
    def from_reference_sequence(cls, ref_seq: ReferenceSequence, pam_variants_: Iterable[PamVariant]) -> PamProtectedReferenceSequence:
        pam_variants: List[PamVariant] = sorted(
            pam_variants_, key=lambda x: x.genomic_position.position)
        pam_seq: str = apply_pam_variants(ref_seq, pam_variants)

        return cls(
            ref_seq.sequence,
            ref_seq.genomic_range,
            pam_seq,
            pam_variants)

    def get_subsequence(self, genomic_range: GenomicRange) -> PamProtectedReferenceSequence:
        ref_seq: ReferenceSequence = super().get_subsequence(genomic_range)
        pam_variants = [
            variant
            for variant in self.pam_variants
            if genomic_range.contains_position(variant.genomic_position)
        ]
        return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_variants)

    def apply_variant(self, variant: BaseVariant, ref_check: bool = False) -> str:
        if not self.genomic_range.contains_position(variant.genomic_position):
            raise ValueError("Variant not in genomic range!")
        # offset: int = variant.genomic_position.position - self.genomic_range.start
        return variant.mutate(self.pam_protected_sequence, self.genomic_range.start, ref_check=ref_check)


@dataclass(frozen=True)
class PamVariant(SubstitutionVariant):
    __slots__ = {'genomic_position', 'ref', 'alt'}

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.ref) != 1 or len(self.alt) != 1:
            raise ValueError(
                "Only single-nucleotide substitutions are allowed "
                "for the purposes of PAM protection!")

    @classmethod
    def from_substitution(cls, var: SubstitutionVariant) -> PamVariant:
        return PamVariant(var.genomic_position, var.ref, var.alt)

    @classmethod
    def from_variant_record(cls, r: VariantRecord) -> PamVariant:
        var: BaseVariant = get_variant(r)
        if not isinstance(var, SubstitutionVariant):
            raise ValueError("PAM protection variants must be substitutions!")
        return cls.from_substitution(var)


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
        return frozenset(self._variants[sgrna_id])

    def get_sgrna_variants_bulk(self, sgrna_ids: FrozenSet[str]) -> FrozenSet[PamVariant]:
        try:
            return frozenset(chain.from_iterable(
                self.get_sgrna_variants(sgrna_id)
                for sgrna_id in sgrna_ids
            ))

        except KeyError as ex:
            sgrna_id: str = ex.args[0]
            raise RuntimeError(f"sgRNA ID '{sgrna_id} not found!'")

    def load(self, fp: str) -> None:

        # Load variants from VCF
        with get_vcf(fp) as variant_file:
            for record in variant_file.fetch():
                sgrna_id: str = record.info['SGRNA'].strip()
                if sgrna_id in self._variants:
                    self._variants[sgrna_id].add(
                        PamVariant.from_variant_record(record))

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
    pam_variants: Set[PamVariant]
) -> PamProtectedReferenceSequence:
    return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_variants)
