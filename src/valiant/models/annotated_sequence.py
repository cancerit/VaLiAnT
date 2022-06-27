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

import abc
from dataclasses import dataclass
from itertools import groupby
import logging
from typing import Dict, List, Tuple

from valiant.models.base import GenomicPosition, StrandedPositionRange
from valiant.enums import MutationType
from valiant.models.codon_table import CodonTable
from valiant.models.variant import SubstitutionVariant


@dataclass(frozen=True)
class BaseAnnotatedSequencePair(abc.ABC):
    pos_range: StrandedPositionRange
    ref_seq: str

    @property
    @abc.abstractmethod
    def alt_seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def ext_ref_seq(self) -> str:
        return self.ref_seq

    @property
    @abc.abstractmethod
    def ext_alt_seq(self) -> str:
        return self.alt_seq

    @property
    @abc.abstractmethod
    def variants(self) -> List[SubstitutionVariant]:
        pass

    @property
    def _ext_offset(self) -> int:
        """Offset converting genomic positions into extended sequence indices"""
        return -self.pos_range.start

    @property
    def variant_positions(self) -> List[GenomicPosition]:
        return [x.genomic_position for x in self.variants]

    def _get_codon_indices(self, positions: List[GenomicPosition]) -> List[int]:
        offset = self._ext_offset
        return [(x.position + offset) // 3 for x in positions]

    def get_codon_indices(self, positions: List[GenomicPosition], **kwargs) -> List[int]:
        return self._get_codon_indices(positions)

    def get_codon_mutation_types_at_codons(self, codon_table: CodonTable, codon_indices: List[int]) -> List[MutationType]:
        return codon_table.get_mutation_types_at(
            self.pos_range.strand,
            self.ext_ref_seq,
            self.ext_alt_seq,
            codon_indices)

    def get_codon_mutation_types_at(self, codon_table: CodonTable, positions: List[GenomicPosition], **kwargs) -> List[MutationType]:
        return self.get_codon_mutation_types_at_codons(
            codon_table, self.get_codon_indices(positions, **kwargs))

    def get_variant_mutation_types(self, codon_table: CodonTable, **kwargs) -> List[MutationType]:
        return self.get_codon_mutation_types_at(
            codon_table, self.variant_positions, **kwargs)


@dataclass(frozen=True)
class ReferenceAnnotatedSequencePair(BaseAnnotatedSequencePair):
    """Non-mutated sequence"""

    __slots__ = ['pos_range', 'ref_seq']

    @property
    def alt_seq(self) -> str:
        return self.ref_seq

    @property
    def ext_ref_seq(self) -> str:
        return self.ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.alt_seq

    @property
    def variants(self) -> List[SubstitutionVariant]:
        return []


@dataclass(frozen=True)
class AnnotatedSequencePair(BaseAnnotatedSequencePair):
    """Potentially mutated sequence (and its reference and variants)"""

    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants']

    _alt_seq: str
    _variants: List[SubstitutionVariant]

    def __post_init__(self) -> None:
        if len(self.ref_seq) != len(self.alt_seq):
            raise ValueError("Mismatching paired sequence lengths!")

    def variant_count(self) -> int:
        return len(self.variants)

    @property
    def variants(self) -> List[SubstitutionVariant]:
        return self._variants

    @property
    def alt_seq(self) -> str:
        return self._alt_seq

    @property
    def ext_ref_seq(self) -> str:
        return self.ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.alt_seq


@dataclass(frozen=True)
class CDSAnnotatedSequencePair(AnnotatedSequencePair):
    """Potentially mutated CDS sequence (and its reference and variants)"""

    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants', 'cds_prefix', 'cds_suffix']

    cds_prefix: str
    cds_suffix: str

    @property
    def cds_ref_seq(self) -> str:
        return f"{self.cds_prefix}{self.ref_seq}{self.cds_suffix}"

    @property
    def cds_alt_seq(self) -> str:
        return f"{self.cds_prefix}{self.alt_seq}{self.cds_suffix}"

    @property
    def ext_ref_seq(self) -> str:
        return self.cds_ref_seq

    @property
    def ext_alt_seq(self) -> str:
        return self.cds_alt_seq

    @property
    def frame(self) -> int:
        return len(self.cds_prefix)

    @property
    def _ext_offset(self) -> int:
        """Offset converting genomic positions into extended sequence indices"""

        return self.frame + super()._ext_offset

    def _log_same_codon_variants(self) -> None:
        offset = self._ext_offset

        def get_variant_codon(variant: SubstitutionVariant) -> int:
            return (variant.genomic_position.position - offset) // 3

        def sort_key(t: Tuple[int, SubstitutionVariant]) -> int:
            return t[0]

        for codon_index, variants in groupby(sorted([
            (get_variant_codon(variant), variant)
            for variant in self.variants
        ], key=sort_key), key=sort_key):
            if len(list(variants)) > 1:
                logging.error("Variants at %s affect the same codon (%d)!" % (
                    ', '.join([str(x.genomic_position) for _, x in variants]),
                    codon_index + 1))

    def get_codon_indices(self, positions: List[GenomicPosition], **kwargs) -> List[int]:
        codon_indices = self._get_codon_indices(positions)
        no_duplicate_codons = kwargs.get('no_duplicate_codons', False)

        # Verify no two variants affect the same codon
        if no_duplicate_codons and len(set(codon_indices)) != len(codon_indices):
            self._log_same_codon_variants()
            raise ValueError("Variant affecting the same codon!")

        return codon_indices
