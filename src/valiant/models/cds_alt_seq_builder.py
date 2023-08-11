########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

from dataclasses import dataclass
from itertools import groupby
import logging
from typing import Dict, List, Set, Tuple

from .alt_seq_builder import AltSeqBuilder
from .base import GenomicPosition, GenomicRange, StrandedPositionRange
from .codon_table import CODON_LENGTH, CodonTable
from .dna_str import DnaStr
from .variant import BaseVariantT
from .variant_group import VariantGroup
from ..enums import MutationType
from ..sgrna_utils import get_codon_index, get_codon_indices_in_range
from ..utils import does_any_set_intersect, has_duplicates


@dataclass(frozen=True)
class CdsAltSeqBuilder(AltSeqBuilder):
    __slots__ = ['gr', 'ref_seq', 'variant_groups', 'cds_prefix', 'cds_suffix']

    cds_prefix: DnaStr
    cds_suffix: DnaStr

    def _extend(self, s: str) -> str:
        return f"{self.cds_prefix}{s}{self.cds_suffix}"

    @property
    def ext_ref_seq(self) -> str:
        return self._extend(self.ref_seq)

    def __post_init__(self) -> None:
        if (
            self.cds_prefix_length > 2 or
            self.cds_suffix_length > 2 or
            self.ext_seq_length % CODON_LENGTH != 0
        ):
            raise ValueError("Invalid extended CDS sequence length: partial codons!")

    @property
    def seq_length(self) -> int:
        return len(self.ref_seq)

    @property
    def ext_seq_length(self) -> int:
        return self.cds_prefix_length + self.cds_suffix_length + self.seq_length

    @property
    def cds_prefix_length(self) -> int:
        return len(self.cds_prefix)

    @property
    def cds_suffix_length(self) -> int:
        return len(self.cds_suffix)

    @property
    def frame(self) -> int:
        return self.cds_prefix_length

    @property
    def codon_count(self) -> int:
        return self.ext_seq_length // CODON_LENGTH

    @property
    def last_codon_index(self) -> int:
        return self.codon_count - 1

    @classmethod
    def from_builder(cls, b: AltSeqBuilder, prefix: str, suffix: str) -> 'CdsAltSeqBuilder':
        return CdsAltSeqBuilder(
            b.gr,
            b.ref_seq,
            b.variant_groups,
            DnaStr(prefix),
            DnaStr(suffix))

    def log_same_codon_variants(self, variant_group_index: int) -> None:

        def get_variant_codon(variant: BaseVariantT) -> int:
            return self.get_ref_codon_index(variant.genomic_position.position)

        def sort_key(t: Tuple[int, BaseVariantT]) -> int:
            return t[0]

        for codon_index, variants in groupby(sorted([
            (get_variant_codon(variant), variant)
            for variant in self.get_variants(variant_group_index)
        ], key=sort_key), key=sort_key):
            if len(list(variants)) > 1:
                logging.warning("Variants at %s affect the same codon (%d)!" % (
                    ', '.join([str(x.genomic_position) for _, x in variants]),
                    codon_index + 1))

    def get_ref_codon_index(self, position: int) -> int:
        """Convert a genomic position into a codon index"""

        return get_codon_index(self.frame, self.start, position)

    def get_ref_codon_indices(self, positions: List[GenomicPosition], no_duplicate_codons: bool = False) -> List[int]:
        # TODO: reintroduce `no_duplicate_codons` validation
        return [self.get_ref_codon_index(x.position) for x in positions]

    def contains_same_codon_variants(self, variant_group_index: int) -> bool:
        codon_indices = self.get_ref_codon_indices(
            self.get_variant_positions(variant_group_index))
        return has_duplicates(codon_indices)

    def get_variant_group_codon_indices(self, variant_group_index: int, **kwargs) -> List[int]:
        positions = self.get_variant_positions(variant_group_index)
        codon_indices = self.get_ref_codon_indices(positions)
        no_duplicate_codons = kwargs.get('no_duplicate_codons', False)

        # Verify no two variants affect the same codon
        if no_duplicate_codons and has_duplicates(codon_indices):
            self.log_same_codon_variants(variant_group_index)
            raise ValueError("Variant affecting the same codon!")

        return codon_indices

    def get_codon_indices_in_range(self, spr: StrandedPositionRange) -> List[int]:
        """Get the indices of the codons spanned by the input range, if any"""

        return get_codon_indices_in_range(
            self.frame,
            self.gr.start,
            self.gr.end,
            spr.start,
            spr.end)

    def get_codon_mutation_types_at_codons(self, codon_table: CodonTable, codon_indices: List[int]) -> List[MutationType]:
        return codon_table.get_mutation_types_at(
            self.gr.strand,
            self.ext_ref_seq,
            self.ext_alt_seq,
            codon_indices)

    def get_codon_mutation_types_at(self, codon_table: CodonTable, positions: List[GenomicPosition], **kwargs) -> List[MutationType]:
        return self.get_codon_mutation_types_at_codons(
            codon_table, self.get_ref_codon_indices(positions, **kwargs))

    def get_variant_mutation_types(self, codon_table: CodonTable, variant_group_index: int, **kwargs) -> List[MutationType]:
        return self.get_codon_mutation_types_at(
            codon_table, self.get_variant_positions(variant_group_index), **kwargs)

    def _get_alt_codon(self, i: int) -> str:
        """Get the ALT sequence of the i-th codon, including the CDS extensions"""

        return self.ext_alt_seq[i * CODON_LENGTH:i * CODON_LENGTH + CODON_LENGTH]

    def get_indexed_alt_codons(self) -> Dict[int, str]:
        """Map codon indices to ALT codons, excluding the CDS extensions"""

        # Map codon indices to ALT codons
        d: Dict[int, str] = {
            i: self._get_alt_codon(i)
            for i in range(self.codon_count)
        }

        # Trim the CDS extension from the first codon, if any
        if self.frame > 0:
            d[0] = d[0][self.frame:]

        # Trim the CDS extension from the last codon, if any
        if self.cds_suffix_length > 0:
            d[self.last_codon_index] = d[self.last_codon_index][:-self.cds_suffix_length]

        return d

    def _get_variant_group_codon_indices(self, g: VariantGroup) -> Set[int]:
        return set().union(*[
            self.get_codon_indices_in_range(variant.ref_range)
            for variant in g.variants
        ])

    def variant_group_codon_clash(self, variant_group_indices: List[int]) -> bool:
        self.validate_variant_group_indices(variant_group_indices)

        if self.variant_group_count < 2 or len(variant_group_indices) < 2:
            return False

        return does_any_set_intersect(
            self._get_variant_group_codon_indices(variant_group)
            for variant_group in (
                self.variant_groups[i]
                for i in variant_group_indices
            )
        )

    def is_variant_nonsynonymous(self, codon_table: CodonTable, variant: BaseVariantT) -> bool:

        # Frame-shifting
        if variant.is_frame_shifting:
            return True

        # Get liminal codons
        start_codon_index = self.get_ref_codon_index(variant.start)
        end_codon_index = self.get_ref_codon_index(variant.ref_end)
        codon_indices = list(range(start_codon_index, end_codon_index + 1))

        # Apply variant to reference sequence
        alt_seq = variant.mutate(self.ref_seq, self.start, ref_check=False)

        # Predict variant effect on the CDS
        mut_types = codon_table.get_mutation_types_at(
            self.gr.strand, self.ref_seq, alt_seq, codon_indices)

        return any(m != MutationType.SYNONYMOUS for m in mut_types)

    def get_sub(self, r: GenomicRange) -> 'CdsAltSeqBuilder':
        raise NotImplementedError
