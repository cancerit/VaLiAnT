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
from typing import Generic, List, Sized, Tuple, TypeVar

from valiant.models.base import GenomicPosition, StrandedPositionRange
from valiant.enums import MutationType
from valiant.models.codon_table import CodonTable
from valiant.models.variant import SubstitutionVariant
from valiant.utils import has_duplicates


RangeT = TypeVar('RangeT', bound='StrandedPositionRange')
VariantT = TypeVar('VariantT', bound='SubstitutionVariant')
AnnotatedSequenceT = TypeVar('AnnotatedSequenceT', bound='BaseAnnotatedSequencePair')


@dataclass(frozen=True)
class BaseAnnotatedSequencePair(abc.ABC, Sized, Generic[RangeT, VariantT]):
    pos_range: RangeT
    ref_seq: str

    def __len__(self) -> int:
        return len(self.ref_seq)

    @property
    @abc.abstractmethod
    def alt_seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def ext_ref_seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def ext_alt_seq(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def variants(self) -> List[VariantT]:
        pass

    @property
    @abc.abstractmethod
    def variant_count(self) -> int:
        pass

    @property
    def variant_positions(self) -> List[GenomicPosition]:
        return [x.genomic_position for x in self.variants]

    @property
    def _ext_offset(self) -> int:
        """Offset converting genomic positions into extended sequence indices"""
        return -self.pos_range.start

    def _get_codon_index(self, position: int) -> int:
        """Convert a genomic position into a codon index"""

        return (position + self._ext_offset) // 3

    def _get_codon_indices(self, positions: List[GenomicPosition]) -> List[int]:
        return [self._get_codon_index(x.position) for x in positions]

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

    def get_variant_codon_indices(self, **kwargs) -> List[int]:
        return self.get_codon_indices(self.variant_positions, **kwargs)

    def get_variant_mutation_types(self, codon_table: CodonTable, **kwargs) -> List[MutationType]:
        return self.get_codon_mutation_types_at(
            codon_table, self.variant_positions, **kwargs)


@dataclass(frozen=True)
class AnnotatedSequencePair(BaseAnnotatedSequencePair, Generic[VariantT]):
    """Potentially mutated sequence (and its reference and variants)"""

    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants']

    _alt_seq: str
    _variants: List[VariantT]

    def __post_init__(self) -> None:
        n: int = len(self.ref_seq)
        if n != len(self.pos_range):
            raise ValueError("Mismatching position range and sequence length!")
        if n != len(self.alt_seq):
            raise ValueError("Mismatching paired sequence lengths!")

    @property
    def variant_count(self) -> int:
        return len(self.variants)

    @property
    def variants(self) -> List[VariantT]:
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
class CDSAnnotatedSequencePair(AnnotatedSequencePair, Generic[VariantT]):
    """Potentially mutated CDS sequence (and its reference and variants)"""

    __slots__ = ['pos_range', 'ref_seq', '_alt_seq', '_variants', 'cds_prefix', 'cds_suffix']

    cds_prefix: str
    cds_suffix: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.ext_seq_length % 3:
            raise ValueError("Invalid length for in-frame sequence!")

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

    @property
    def ext_seq_length(self) -> int:
        return len(self) + self.cds_prefix_length + self.cds_suffix_length

    @property
    def codon_count(self) -> int:
        return self.ext_seq_length // 3

    @property
    def last_codon_index(self) -> int:
        return self.codon_count - 1

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
    def start(self) -> int:
        return self.pos_range.start

    @property
    def end(self) -> int:
        return self.pos_range.end

    @property
    def cds_sequence_start(self) -> int:
        return self.start - self.frame

    @property
    def cds_sequence_end(self) -> int:
        return self.end + self.cds_suffix_length

    @property
    def _ext_offset(self) -> int:
        """Offset converting genomic positions into extended sequence indices"""

        return self.frame + super()._ext_offset

    def log_same_codon_variants(self) -> None:
        offset = self._ext_offset

        def get_variant_codon(variant: VariantT) -> int:
            return (variant.genomic_position.position - offset) // 3

        def sort_key(t: Tuple[int, VariantT]) -> int:
            return t[0]

        for codon_index, variants in groupby(sorted([
            (get_variant_codon(variant), variant)
            for variant in self.variants
        ], key=sort_key), key=sort_key):
            if len(list(variants)) > 1:
                logging.warning("Variants at %s affect the same codon (%d)!" % (
                    ', '.join([str(x.genomic_position) for _, x in variants]),
                    codon_index + 1))

    @property
    def contains_same_codon_variants(self) -> bool:
        codon_indices = super().get_codon_indices(self.variant_positions)
        return has_duplicates(codon_indices)

    def get_codon_indices(self, positions: List[GenomicPosition], **kwargs) -> List[int]:
        codon_indices = super().get_codon_indices(positions)
        no_duplicate_codons = kwargs.get('no_duplicate_codons', False)

        # Verify no two variants affect the same codon
        if no_duplicate_codons and has_duplicates(codon_indices):
            self.log_same_codon_variants()
            raise ValueError("Variant affecting the same codon!")

        return codon_indices

    def get_codon_indices_in_range(self, spr: StrandedPositionRange) -> List[int]:
        """Get the indices of the codons spanned by the input range, if any"""

        # Check an intersection exists
        # NOTE: `pos_range` may be a GenomicRange, in which case __contains__
        # would fail due to the absence of the chromosome attribute in `spr`

        # Fail on incorrect strand (pathological state)
        if spr.strand != self.pos_range.strand:
            raise ValueError("Failed to retrieve codon indices: incorrect strand!")

        start: int = self.pos_range.start
        end: int = self.pos_range.end

        # NOTE: in a `PositionRange`, `end` is guaranteed to be larger than `start`
        if spr.end < start or spr.start > end:
            return []

        # Get first codon index
        codon_start: int = 0 if spr.start <= start else self._get_codon_index(spr.start)  # (spr.start + self._ext_offset) // 3  # (spr.start - self.cds_sequence_start) // 3

        if spr.end == spr.start:
            return [codon_start]

        # Get last codon index
        codon_end: int = (
            self.last_codon_index if spr.end >= end else
            self._get_codon_index(spr.end)
        )

        # Generate full codon index range
        return (
            [codon_start] if codon_end == codon_start else
            list(range(codon_start, codon_end + 1))
        )
