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
import csv
from dataclasses import dataclass
from functools import lru_cache
from itertools import chain
import logging
import re
from typing import Dict, Iterable, List, Optional, FrozenSet, Tuple
import pandas as pd
from pyranges import PyRanges
from .base import GenomicRange
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


@dataclass(init=False)
class ReferenceSequenceRanges:
    __slots__ = {
        'ref_range',
        'sgrna_ids',
        '_const_regions',
        '_target_regions'
    }

    ref_range: GenomicRange
    sgrna_ids: FrozenSet[str]
    _const_regions: Tuple[Optional[GenomicRange], Optional[GenomicRange]]
    _target_regions: Tuple[
        Optional[TargetReferenceRegion],
        TargetReferenceRegion,
        Optional[TargetReferenceRegion]
    ]

    def __init__(
        self,
        chromosome: str,
        strand: str,
        ref_start: int,
        ref_end: int,
        target_region_2_start: int,
        target_region_2_end: int,
        target_region_2_extension: Tuple[int, int],
        mutators: Tuple[FrozenSet[TargetonMutator], FrozenSet[TargetonMutator], FrozenSet[TargetonMutator]],
        sgrna_ids: FrozenSet[str]
    ) -> None:

        def get_genomic_range(start: int, end: int) -> GenomicRange:
            return GenomicRange(chromosome, start, end, strand)

        self.sgrna_ids = sgrna_ids

        if len(target_region_2_extension) != 2 or any(x < 0 for x in target_region_2_extension):
            raise ValueError("Invalid extension vector!")

        # Get lenghts of target regions 1 and 3
        r1_len, r3_len = target_region_2_extension

        self.ref_range = get_genomic_range(ref_start, ref_end)
        target_region_2_range = get_genomic_range(
            target_region_2_start, target_region_2_end)

        if target_region_2_range not in self.ref_range:
            raise ValueError("Target region 2 is outside of the reference genomic range!")

        if (
            target_region_2_range.start - r1_len < self.ref_range.start
            or target_region_2_range.end + r3_len > self.ref_range.end
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
            self.ref_range.start,
            targeton_start - 1
        ) if targeton_start - self.ref_range.start > 0 else None

        const_region_2: Optional[GenomicRange] = get_genomic_range(
            targeton_end + 1,
            self.ref_range.end
        ) if self.ref_range.end - targeton_end > 0 else None

        # Store regions
        self._const_regions = (const_region_1, const_region_2)
        self._target_regions = (
            TargetReferenceRegion(target_region_1_range, mutators[0]) if target_region_1_range else None,
            TargetReferenceRegion(target_region_2_range, mutators[1]),
            TargetReferenceRegion(target_region_3_range, mutators[2]) if target_region_3_range else None
        )

    @staticmethod
    def parse_mutator_tuples(s: str) -> List[FrozenSet[TargetonMutator]]:
        m: Optional[re.Match] = mutator_vector_re.match(s)

        if not m:
            raise ValueError("Invalid format for vector!")

        return [
            parse_mutators(mutator_group) if mutator_group else set()
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
        sgrna_ids: FrozenSet[str] = set(parse_list(row[8]))

        return cls(
            row[0],
            row[1],
            int(row[2]),
            int(row[3]),
            int(row[4]),
            int(row[5]),
            (extensions[0], extensions[1]),
            (mutators[0], mutators[1], mutators[2]),
            sgrna_ids)

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
    __slots__ = {'_rsrs', '_ref_ranges', '_region_ranges'}

    _rsrs: Dict[int, ReferenceSequenceRanges]
    _ref_ranges: PyRanges
    _region_ranges: PyRanges

    def __init__(self, rsrs: Iterable[ReferenceSequenceRanges]) -> None:

        # Assign numerical indices to targetons
        self._rsrs = {rsr_id: rsr for rsr_id, rsr in enumerate(rsrs)}

        # Collect reference sequence unstranded genomic ranges (unique)
        self._ref_ranges = genomic_ranges_to_unstranded_pyranges(
            rsr.ref_range for rsr in self._rsrs.values())

        # Collect subregion unstranded genomic ranges
        df: pd.DataFrame = pd.DataFrame.from_records(chain([
            (*r, i)
            for i, rsr in self._rsrs.items()
            for r in rsr.regions_as_pyranges()
        ]), columns=['Chromosome', 'Strand', 'Start', 'End', 'is_const', 'targeton_id'])
        df.targeton_id = df.targeton_id.astype(get_smallest_int_type(max(self._rsrs.keys())))
        self._region_ranges = PyRanges(df=df)

    @classmethod
    def load(cls, fp: str) -> ReferenceSequenceRangeCollection:
        return cls(map(ReferenceSequenceRanges.from_row, load_tsv(fp, CSV_HEADER)))

    @property
    def sgrna_ids(self) -> FrozenSet[str]:
        return frozenset.union(*[rsr.sgrna_ids for rsr in self._rsrs.values()])

    @property
    def target_ranges(self) -> PyRanges:
        return self._region_ranges[~self._region_ranges.is_const]

    @property
    def ref_ranges(self) -> FrozenSet[GenomicRange]:
        return frozenset(rsr.ref_range for rsr in self._rsrs.values())

    @property
    def strands(self) -> FrozenSet[str]:
        return frozenset(gr.strand for gr in self.ref_ranges)

    @property
    def mutarors(self) -> FrozenSet[TargetonMutator]:
        return frozenset.union(*[rsr.mutators for rsr in self._rsrs.values()])
