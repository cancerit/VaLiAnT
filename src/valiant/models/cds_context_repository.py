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
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from pyranges import PyRanges

from .base import GenomicRange, GenomicRangePair, TranscriptInfo
from .cds_info import CdsInfo
from .exon_info import ExonInfo
from .exon_repository import ExonRepository


def compute_cds_contexts(targets: PyRanges, exons: ExonRepository) -> Dict[GenomicRange, Tuple[ExonInfo, GenomicRangePair]]:

    def get_delta_start(df: pd.DataFrame, mask: pd.Series) -> pd.Series:
        return df.loc[mask, 'Start'] - df.loc[mask, 'Start_b']

    def get_delta_end(df: pd.DataFrame, mask: pd.Series) -> pd.Series:
        return df.loc[mask, 'End_b'] - df.loc[mask, 'End']

    exonic_ranges: pd.DataFrame = targets.join(
        exons.exon_ranges, strandedness='same').drop(['Strand_b']).as_df()

    # Skip if no CDS targets are found
    if exonic_ranges.shape[0] == 0:
        return {}

    # Compute target offsets within exons
    plus_strand_mask: pd.Series = exonic_ranges.Strand == '+'
    minus_strand_mask: pd.Series = ~plus_strand_mask
    exonic_ranges.loc[plus_strand_mask, 'delta_5p'] = get_delta_start(exonic_ranges, plus_strand_mask)
    exonic_ranges.loc[minus_strand_mask, 'delta_5p'] = get_delta_end(exonic_ranges, minus_strand_mask)
    exonic_ranges.loc[plus_strand_mask, 'delta_3p'] = get_delta_end(exonic_ranges, plus_strand_mask)
    exonic_ranges.loc[minus_strand_mask, 'delta_3p'] = get_delta_start(exonic_ranges, minus_strand_mask)
    del plus_strand_mask, minus_strand_mask

    exonic_ranges['exon_frame'] = exonic_ranges.frame
    exonic_ranges['len'] = exonic_ranges.End - exonic_ranges.Start
    exonic_ranges.frame = (exonic_ranges.exon_frame + exonic_ranges.delta_5p.mod(3).astype(np.int8)).mod(3)

    # Detect and validate partial exon targets
    if exonic_ranges.delta_5p.lt(0).any() or exonic_ranges.delta_3p.lt(0).any():
        raise NotImplementedError("Mixed exonic and intronic targets not supported!")

    exonic_ranges['cds_ext_3_length'] = (3 - (exonic_ranges.len + exonic_ranges.frame) % 3) % 3

    target_cds_extension_info = [
        CdsInfo.from_pyr(k, record)
        for k, ranges in PyRanges(df=exonic_ranges)
        for record in ranges.itertuples()
    ]
    del exonic_ranges

    # TODO: optimise the retrieval of exons by transcript and index
    return {
        x.genomic_range: (
            x.get_exon_info(),
            exons.get_exon_ext_genomic_ranges(x))
        for x in target_cds_extension_info
    }


@dataclass(frozen=True)
class CDSContextRepository:
    """
    Positional and reading frame information on the transcript exons

    - cds_ranges: genomic ranges of all transcript exons
    - _target_ranges: regions of interest (to filter the exons)
    - _target_cds_contexts: non-adjacent nucleotide genomic ranges of the
      liminal codons of the exons of interest
    """

    __slots__ = ['target_ranges', '_target_cds_contexts']

    target_ranges: PyRanges
    _target_cds_contexts: Dict[GenomicRange, Tuple[ExonInfo, GenomicRangePair]]

    @classmethod
    def from_exons(cls, targets: PyRanges, exons: ExonRepository) -> CDSContextRepository:
        return cls(targets, compute_cds_contexts(targets, exons))

    def get_cds_extensions(self, genomic_range: GenomicRange) -> Optional[GenomicRangePair]:
        if not self._target_cds_contexts or genomic_range not in self._target_cds_contexts:
            return None
        return self._target_cds_contexts[genomic_range][1]

    def get_all_cds_extensions(self) -> Set[GenomicRange]:
        return set(
            genomic_range
            for genomic_range in chain.from_iterable((
                cds_ranges
                for exon_info, cds_ranges in
                self._target_cds_contexts.values()))
            if genomic_range
        ) if self._target_cds_contexts else set()

    def get_exon_info(self, genomic_range: GenomicRange) -> Optional[ExonInfo]:
        if not self._target_cds_contexts or genomic_range not in self._target_cds_contexts:
            return None
        return self._target_cds_contexts[genomic_range][0]

    def get_transcript_info(self, genomic_range: GenomicRange) -> Optional[TranscriptInfo]:
        exon_info: Optional[ExonInfo] = self.get_exon_info(genomic_range)
        return exon_info.transcript_info if exon_info else None

    def get_exons(self) -> List[ExonInfo]:
        return [
            exon_info
            for exon_info, _ in self._target_cds_contexts.values()
        ]
