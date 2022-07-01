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
import os
from typing import Any
import pandas as pd
from ..constants import METADATA_FIELDS_SET, METADATA_FIELDS
from ..loaders.vcf import write_vcf
from ..utils import get_constant_category
from ..writers import write_oligo_metadata, write_oligo_unique
from .oligo_generation_info import OligoGenerationInfo
from .refseq_repository import ReferenceSequenceRepository
from .variant import get_records


def _fill_metadata(species: str, assembly: str, metadata: pd.DataFrame) -> pd.DataFrame:

    # Add global metadata
    rown: int = metadata.shape[0]
    metadata['species'] = get_constant_category(species, rown)
    metadata['assembly'] = get_constant_category(assembly, rown)

    # Add missing columns
    for field in METADATA_FIELDS_SET - set(metadata.columns):
        metadata[field] = None

    return metadata


@dataclass(frozen=True)
class MetadataTable:
    __slots__ = {'metadata', '_oligo_length_mask', '_oligo_min_length_mask', '_oligo_max_length_mask', 'too_short_oligo_n', 'short_oligo_n', 'long_oligo_n'}

    metadata: pd.DataFrame
    _oligo_length_mask: pd.Series
    _oligo_min_length_mask: pd.Series
    _oligo_max_length_mask: pd.Series
    too_short_oligo_n: int
    short_oligo_n: int
    long_oligo_n: int

    @property
    def oligo_length_mask(self) -> pd.Series:
        return self._oligo_length_mask

    def _setattr(self, attr: str, value: Any) -> None:
        object.__setattr__(self, attr, value)

    def __init__(self, metadata: pd.DataFrame, oligo_min_length: int, oligo_max_length: int) -> None:
        if oligo_max_length < 1:
            raise ValueError("Invalid maximum oligonucleotide length!")

        self._setattr('metadata', metadata)

        # Create olignonucleotide length masks
        self._setattr('_oligo_min_length_mask', metadata.oligo_length >= oligo_min_length)
        self._setattr('_oligo_max_length_mask', metadata.oligo_length <= oligo_max_length)
        self._setattr('_oligo_length_mask', self._oligo_min_length_mask & self._oligo_max_length_mask)

        # Count the oligonucleotides per length range
        n: int = len(self.oligo_length_mask)
        self._setattr('too_short_oligo_n', n - self._oligo_min_length_mask.sum())
        self._setattr('long_oligo_n', n - self._oligo_max_length_mask.sum())
        self._setattr('short_oligo_n', n - self.too_short_oligo_n - self.long_oligo_n)

    @classmethod
    def from_partial(
        cls,
        species: str,
        assembly: str,
        metadata: pd.DataFrame,
        oligo_min_length: int,
        oligo_max_length: int
    ) -> MetadataTable:
        return cls(
            _fill_metadata(species, assembly, metadata),
            oligo_min_length,
            oligo_max_length)

    def write_metadata_file(self, fp: str) -> None:
        write_oligo_metadata((
            self.metadata.loc[self.oligo_length_mask, METADATA_FIELDS] if self.long_oligo_n > 0 else
            self.metadata[METADATA_FIELDS]
        ), fp)

    def write_excluded_metadata_file(self, fp: str) -> None:
        write_oligo_metadata(
            self.metadata.loc[~self.oligo_length_mask, METADATA_FIELDS], fp)

    def write_unique_file(self, fp: str) -> None:
        write_oligo_unique(
            self.metadata.loc[
                self.oligo_length_mask,
                ['oligo_name', 'mseq']
            ].sort_values(
                'oligo_name',
                kind='mergesort'
            ).drop_duplicates(subset=['mseq'], keep='first'),
            fp)

    def write_vcf_file(self, fp: str, ref_repository: ReferenceSequenceRepository) -> None:
        metadata = (
            self.metadata[self.oligo_length_mask] if self.long_oligo_n > 0 else
            self.metadata
        )
        write_vcf(
            fp,
            list(metadata.ref_chr.cat.categories.values),
            get_records(ref_repository, metadata))

    def write_common_files(self, out_dir: str, base_fn: str) -> None:
        if self.short_oligo_n > 0:

            # Save metadata to file (filtered and reordered columns)
            metadata_fn = base_fn + '_meta.csv'
            self.write_metadata_file(os.path.join(out_dir, metadata_fn))

            # Save unique oligonucleotides to file (and filter by length)
            unique_oligos_fn = base_fn + '_unique.csv'
            self.write_unique_file(os.path.join(out_dir, unique_oligos_fn))

        if self.long_oligo_n > 0:

            # Save discarded metadata to file
            excluded_metadata_fn = base_fn + '_meta_excluded.csv'
            self.write_excluded_metadata_file(
                os.path.join(out_dir, excluded_metadata_fn))

    def write_sge_files(
        self,
        out_dir: str,
        base_fn: str,
        ref_repository: ReferenceSequenceRepository
    ) -> None:
        self.write_common_files(out_dir, base_fn)
        if self.short_oligo_n > 0:

            # Save variants to file (VCF format)
            vcf_fn = base_fn + '.vcf'
            self.write_vcf_file(os.path.join(out_dir, vcf_fn), ref_repository)

    def get_info(self) -> OligoGenerationInfo:
        return OligoGenerationInfo(self.long_oligo_n)
