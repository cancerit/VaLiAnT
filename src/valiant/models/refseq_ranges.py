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
import re
from typing import Dict, Iterable, List, Optional, FrozenSet, Tuple
import pandas as pd
from pyranges import PyRanges
from .base import GenomicRange, PositionRange
from ..constants import PYR_CHR, PYR_END, PYR_START, PYR_STRAND
from ..enums import TargetonMutator
from ..loaders.tsv import load_tsv
from ..utils import get_smallest_int_type, parse_list, parse_mutators


CSV_HEADER = [
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'r2_start',
    'r2_end',
    'ext_vector',
    'action_vector',
    'sgrna_vector'
]

# Mutation vector pattern, e.g.: `(1del), (snv, 1del), (3del)`
mutator_vector_re: re.Pattern = re.compile(
    r'\s*,\s*'.join([r'\((\s*[\w\-_]+\s*(?:,\s*[\w\-_]+)*)?\s*\)'] * 3))


def genomic_ranges_to_pyranges(genomic_ranges: Iterable[GenomicRange]) -> PyRanges:
    return PyRanges(df=pd.DataFrame.from_records(chain([
        genomic_range.as_pyrange()
        for genomic_range in genomic_ranges
    ]), columns=['Chromosome', 'Strand', 'Start', 'End']).drop_duplicates())


def genomic_ranges_to_unstranded_pyranges(genomic_ranges: Iterable[GenomicRange]) -> PyRanges:
    return PyRanges(df=pd.DataFrame.from_records(chain([
        genomic_range.as_unstranded_pyrange()
        for genomic_range in genomic_ranges
    ]), columns=['Chromosome', 'Start', 'End']).drop_duplicates())


@dataclass(frozen=True)
class TargetReferenceRegion:
    __slots__ = {'genomic_range', 'mutators'}

    genomic_range: GenomicRange
    mutators: FrozenSet[TargetonMutator]


@dataclass
class TargetonConfig:
    __slots__ = [
        'chromosome',
        'strand',
        'ref_start',
        'ref_end',
        'target_region_2_start',
        'target_region_2_end',
        'target_region_2_extension',
        'mutators',
        'sgrna_ids'
    ]

    chromosome: str
    strand: str
    ref_start: int
    ref_end: int
    target_region_2_start: int
    target_region_2_end: int
    target_region_2_extension: Tuple[int, int]
    mutators: Tuple[FrozenSet[TargetonMutator], FrozenSet[TargetonMutator], FrozenSet[TargetonMutator]]
    sgrna_ids: FrozenSet[str]


@dataclass
class ReferenceSequenceRanges:
    __slots__ = {
        'ref_range',
        'sgrna_ids',
        '_bg_mask_regions',
        '_const_regions',
        '_target_regions'
    }

    ref_range: GenomicRange
    sgrna_ids: FrozenSet[str]
    _bg_mask_regions: List[GenomicRange]
    _const_regions: Tuple[Optional[GenomicRange], Optional[GenomicRange]]
    _target_regions: Tuple[
        Optional[TargetReferenceRegion],
        TargetReferenceRegion,
        Optional[TargetReferenceRegion]
    ]

    @property
    def bg_ranges(self) -> List[GenomicRange]:
        return self.ref_range.diff(self._bg_mask_regions)

    @classmethod
    def from_config(
        cls,
        config: TargetonConfig
    ) -> ReferenceSequenceRanges:

        # TODO: implement
        bg_mask_regions = []

        def get_genomic_range(start: int, end: int) -> GenomicRange:
            return GenomicRange(config.chromosome, start, end, config.strand)

        if len(config.target_region_2_extension) != 2 or any(x < 0 for x in config.target_region_2_extension):
            raise ValueError("Invalid extension vector!")

        # Get lenghts of target regions 1 and 3
        r1_len, r3_len = config.target_region_2_extension

        ref_range = get_genomic_range(config.ref_start, config.ref_end)
        target_region_2_range = get_genomic_range(
            config.target_region_2_start, config.target_region_2_end)

        if target_region_2_range not in ref_range:
            raise ValueError("Target region 2 is outside of the reference genomic range!")

        if (
            target_region_2_range.start - r1_len < ref_range.start
            or target_region_2_range.end + r3_len > ref_range.end
        ):
            raise ValueError("Invalid extension vector: exceeding reference genomic range!")

        # Calculate adjacent target region genomic ranges
        target_region_1_range: Optional[GenomicRange] = get_genomic_range(
            target_region_2_range.start - r1_len,
            target_region_2_range.start - 1) if r1_len > 0 else None

        target_region_3_range: Optional[GenomicRange] = get_genomic_range(
            target_region_2_range.end + 1,
            target_region_2_range.end + r3_len) if r3_len > 0 else None

        # Calculate targeton genomic range
        targeton_start: int = target_region_1_range.start if target_region_1_range else target_region_2_range.start
        targeton_end: int = target_region_3_range.end if target_region_3_range else target_region_2_range.end

        # Calculate flanking constant region genomic ranges
        const_region_1: Optional[GenomicRange] = get_genomic_range(
            ref_range.start,
            targeton_start - 1
        ) if targeton_start - ref_range.start > 0 else None

        const_region_2: Optional[GenomicRange] = get_genomic_range(
            targeton_end + 1,
            ref_range.end
        ) if ref_range.end - targeton_end > 0 else None

        # Store regions
        const_regions = (const_region_1, const_region_2)
        target_regions = (
            TargetReferenceRegion(target_region_1_range, config.mutators[0]) if target_region_1_range else None,
            TargetReferenceRegion(target_region_2_range, config.mutators[1]),
            TargetReferenceRegion(target_region_3_range, config.mutators[2]) if target_region_3_range else None
        )

        return cls(ref_range, config.sgrna_ids, bg_mask_regions, const_regions, target_regions)

    @staticmethod
    def parse_mutator_tuples(s: str) -> List[FrozenSet[TargetonMutator]]:
        m: Optional[re.Match] = mutator_vector_re.match(s)

        if not m:
            raise ValueError("Invalid format for vector!")

        return [
            parse_mutators(mutator_group) if mutator_group else frozenset()
            for mutator_group in m.groups()
        ]

    @classmethod
    def from_row(cls, row: List[str]) -> ReferenceSequenceRanges:

        # Target region extension vector
        extensions: List[int] = list(map(int, parse_list(row[6])))
        if len(extensions) != 2:
            raise ValueError("Invalid extension vector: two values expected!")

        # Action vector
        mutators: List[FrozenSet[TargetonMutator]] = cls.parse_mutator_tuples(row[7])

        # sgRNA ID vector
        sgrna_ids: FrozenSet[str] = frozenset(parse_list(row[8]))

        return cls.from_config(TargetonConfig(
            row[0],
            row[1],
            int(row[2]),
            int(row[3]),
            int(row[4]),
            int(row[5]),
            (extensions[0], extensions[1]),
            (mutators[0], mutators[1], mutators[2]),
            sgrna_ids))

    @property
    def mutators(self) -> FrozenSet[TargetonMutator]:
        return frozenset.union(*[trr.mutators for trr in self._target_regions if trr])

    @property
    def target_regions(self) -> List[TargetReferenceRegion]:
        return [trr for trr in self._target_regions if trr]

    @property
    def const_region_ranges(self) -> List[GenomicRange]:
        return [gr for gr in self._const_regions if gr]

    @property
    def const_region_1(self) -> Optional[GenomicRange]:
        return self._const_regions[0]

    @property
    def const_region_2(self) -> Optional[GenomicRange]:
        return self._const_regions[1]

    @property
    def target_ragion_1(self) -> Optional[GenomicRange]:
        return self._target_regions[0].genomic_range if self._target_regions[0] else None

    @property
    def target_ragion_2(self) -> GenomicRange:
        return self._target_regions[1].genomic_range

    @property
    def target_ragion_3(self) -> Optional[GenomicRange]:
        return self._target_regions[2].genomic_range if self._target_regions[2] else None

    def is_range_in_constant_region(self, gr: PositionRange) -> bool:
        for cr in self._const_regions:
            if cr is not None and not (gr.start > cr.end or gr.end < cr.start):
                return True
        return False

    def regions_as_pyranges(self) -> List[Tuple[str, str, int, int, bool]]:
        return [
            (*gr.as_pyrange(), True)
            for gr in self.const_region_ranges
        ] + [
            (*trr.genomic_range.as_pyrange(), len(trr.mutators) == 0)
            for trr in self.target_regions
        ]


@dataclass(init=False)
class ReferenceSequenceRangeCollection:
    __slots__ = {'_rsrs', '_ref_ranges', '_ref_bg_ranges', '_region_ranges'}

    _rsrs: Dict[int, ReferenceSequenceRanges]
    _unmasked_ref_ranges: PyRanges
    _ref_ranges: PyRanges
    _ref_bg_ranges: PyRanges
    _region_ranges: PyRanges

    def __init__(self, rsrs: Iterable[ReferenceSequenceRanges]) -> None:

        # Assign numerical indices to targetons
        self._rsrs = {rsr_id: rsr for rsr_id, rsr in enumerate(rsrs)}

        # Collect reference sequence unstranded genomic ranges (unique)
        self._ref_ranges = genomic_ranges_to_unstranded_pyranges(
            rsr.ref_range for rsr in self._rsrs.values())

        # Collect reference sequence unstranded genomic ranges
        # for the purposes of background variant filtering
        self._ref_bg_ranges = genomic_ranges_to_unstranded_pyranges(
            chain.from_iterable(
                rsr.bg_ranges for rsr in self._rsrs.values()))

        # Collect subregion unstranded genomic ranges
        df: pd.DataFrame = pd.DataFrame.from_records(chain([
            (*r, i)
            for i, rsr in self._rsrs.items()
            for r in rsr.regions_as_pyranges()
        ]), columns=[PYR_CHR, PYR_STRAND, PYR_START, PYR_END, 'is_const', 'targeton_id'])
        df.targeton_id = df.targeton_id.astype(get_smallest_int_type(max(self._rsrs.keys())))
        self._region_ranges = PyRanges(df=df)

    @classmethod
    def load(cls, fp: str) -> ReferenceSequenceRangeCollection:
        return cls(map(ReferenceSequenceRanges.from_row, load_tsv(fp, CSV_HEADER)))

    @property
    def sgrna_ids(self) -> FrozenSet[str]:
        return frozenset.union(*[rsr.sgrna_ids for rsr in self._rsrs.values()])

    @property
    def region_ranges(self) -> PyRanges:
        return self._region_ranges

    @property
    def target_ranges(self) -> PyRanges:
        return self._region_ranges[~self._region_ranges.is_const]

    @property
    def ref_ranges(self) -> FrozenSet[GenomicRange]:
        return frozenset(rsr.ref_range for rsr in self._rsrs.values())

    @property
    def ref_pyr(self) -> PyRanges:
        return self._ref_ranges

    @property
    def ref_bg_pyr(self) -> PyRanges:
        return self._ref_bg_ranges

    @property
    def strands(self) -> FrozenSet[str]:
        return frozenset(gr.strand for gr in self.ref_ranges)

    @property
    def mutarors(self) -> FrozenSet[TargetonMutator]:
        return frozenset.union(*[rsr.mutators for rsr in self._rsrs.values()])
