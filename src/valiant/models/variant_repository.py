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
from typing import Callable, Dict, Generic, List, Set, Tuple, TypeVar

import pandas as pd
from pyranges import PyRanges

from ..loaders.vcf import load_vcf, load_vcf_manifest, var_type_sub, var_type_del, var_type_ins, var_class_unclass, var_class_mono
from ..utils import get_id_column, get_var_types
from .base import GenomicPosition, GenomicRange
from .custom_variants import CustomVariant
from .variant import BaseVariant, DeletionVariant, InsertionVariant, SubstitutionVariant


VariantRepositoryVariant = TypeVar('VariantRepositoryVariant', bound=CustomVariant)
VariantConstructor = Callable[[pd.Series], BaseVariant]
VariantRepositoryVariantConstructor = Callable[[pd.Series], VariantRepositoryVariant]


VAR_TYPE_CONSTRUCTOR: Dict[int, VariantConstructor] = {
    var_type_sub: lambda t: SubstitutionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.ref, t.alt),
    var_type_del: lambda t: DeletionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.ref),
    var_type_ins: lambda t: InsertionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.alt)
}


def _map_any_variants(
    variants: pd.DataFrame,
    constructor: VariantRepositoryVariantConstructor
) -> Dict[int, VariantRepositoryVariant]:
    return {
        t.variant_id: constructor(t)
        for t in variants.itertuples(index=False)
    }


def _map_variants(variants: pd.DataFrame, var_type: int, is_custom: bool, mask: bool = True) -> Dict[int, VariantRepositoryVariant]:
    constructor = VAR_TYPE_CONSTRUCTOR[var_type]

    def get_custom_variant(t: pd.Series) -> CustomVariant:
        return CustomVariant(
            constructor(t),
            t.vcf_alias,
            t.vcf_var_id)

    def get_base_variant(t: pd.Series) -> CustomVariant:
        return CustomVariant(constructor(t), None, None)

    return _map_any_variants(
        (
            variants[variants.var_type == var_type] if mask else
            variants
        ), (
            get_custom_variant if is_custom else
            get_base_variant
        )
    )


def _regions_to_chromosome_boundaries(regions: PyRanges) -> Dict[str, Tuple[int, int]]:
    return {
        chromosome: (df.Start.min() + 1, df.End.max())
        for chromosome, df in regions.dfs.items()  # type: ignore
    }


@dataclass
class VariantRepository(Generic[VariantRepositoryVariant]):
    _variants: Dict[int, VariantRepositoryVariant]
    _region_variants: Dict[Tuple[str, int, int], Set[int]]

    @classmethod
    def load_vcf(cls, vcf_fp: str, regions: PyRanges) -> VariantRepository:
        # TODO: verify variant filtering criteria being applied at this stage (if any)
        chromosome_boundaries = _regions_to_chromosome_boundaries(regions)
        df = load_vcf(vcf_fp, chromosome_boundaries)
        # TODO: discuss variant identifiers
        df['vcf_alias'] = None
        df['vcf_var_id'] = None
        return cls.from_df(regions, df)

    @classmethod
    def load(cls, manifest_fp: str, regions: PyRanges) -> VariantRepository:

        # Load all permitted variants from multiple VCF files
        chromosome_boundaries = _regions_to_chromosome_boundaries(regions)
        df = load_vcf_manifest(manifest_fp, chromosome_boundaries)
        return cls.from_df(regions, df)

    @classmethod
    def from_df(cls, regions: PyRanges, custom_variants: pd.DataFrame) -> VariantRepository:
        custom_variants_n: int = custom_variants.shape[0]

        if custom_variants_n == 0:
            return cls({}, {})

        custom_variants_n: int = custom_variants.shape[0]

        logging.debug("Collected %d custom variants." % custom_variants_n)

        # Make start positions zero-based
        custom_variants.Start -= 1

        # Assign identifier to all variants
        custom_variants['variant_id'] = get_id_column(custom_variants_n)
        custom_variant_ranges: PyRanges = PyRanges(df=custom_variants)

        # Match regions with variants
        ref_ranges_variants: pd.DataFrame = regions.join(
            custom_variant_ranges,
            suffix='_var'
        )[[
            'Start_var',
            'End_var',
            'variant_id',
            'ref',
            'alt',
            'vcf_alias',
            'vcf_var_id',
            'var_type',
            'var_class'
        ]].as_df()
        del custom_variant_ranges

        variants: pd.DataFrame = ref_ranges_variants.drop([
            'Start',
            'End'
        ], axis=1).drop_duplicates([
            'variant_id'
        ])

        # Log and discard monomorphic variants
        mono_mask: pd.Series = variants.var_class == var_class_mono
        if mono_mask.any():
            for chromosome, start in variants.loc[mono_mask, [
                'Chromosome',
                'Start_var'
            ]].itertuples(index=False, name=None):
                logging.info(f"Monomorphic variant at {chromosome}:{start + 1} (SKIPPED).")
            variants = variants[~mono_mask]
            mono_mask = ref_ranges_variants.var_class != var_class_mono
            ref_ranges_variants = ref_ranges_variants[mono_mask]
        del mono_mask

        # Log unclassified variants
        unclass_mask: pd.Series = variants.var_class == var_class_unclass
        for r in variants[unclass_mask].itertuples(index=False):
            logging.info(f"Unclassified variant at {r.Chromosome}:{r.Start_var + 1}: {r.ref}>{r.alt}.")
        del unclass_mask

        # List all variant types represented in the collection
        var_types: List[int] = get_var_types(variants.var_type)

        var_types_n: int = len(var_types)
        if not var_types_n:
            raise RuntimeError("No variant types in variant table!")

        # Map variant indices to variant objects
        var_types_gt_one: bool = var_types_n > 1
        matching_variants = _map_variants(
            variants, var_types[0], True, mask=var_types_gt_one)
        if var_types_gt_one:
            for var_type in var_types[1:]:
                matching_variants.update(_map_variants(variants, var_type, True))
        del variants, var_types_gt_one

        # Map regions to variant indices
        ref_ranges_variant_ids: Dict[Tuple[str, int, int], Set[int]] = {
            (chromosome, start + 1, end): set(g.variant_id[(
                (g.Start_var >= start)
                & (g.End_var <= end)
            )].unique())
            for (chromosome, start, end), g in ref_ranges_variants.groupby([
                'Chromosome',
                'Start',
                'End'
            ])
        }

        return cls(matching_variants, ref_ranges_variant_ids)

    def get_variants(self, genomic_range: GenomicRange) -> Set[VariantRepositoryVariant]:
        r: Tuple[str, int, int] = genomic_range.as_unstranded()

        if r not in self._region_variants or not self._region_variants[r]:
            return set()

        return set(self._variants[var_id] for var_id in self._region_variants[r])
