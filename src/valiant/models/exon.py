########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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
from typing import Dict, Optional, Set, Tuple
import numpy as np
import pandas as pd
from pyranges import PyRanges
from .base import GenomicRange, TranscriptInfo
from ..loaders.gff import load_gff_cds
from ..utils import validate_strand

GenomicRangePair = Tuple[Optional[GenomicRange], Optional[GenomicRange]]


def pyrange_to_genomic_range(k: Tuple[str, str], record) -> GenomicRange:
    return GenomicRange(k[0], record.Start + 1, record.End, k[1])


@dataclass
class ExonInfo:
    __slots__ = {'transcript_info', 'genomic_range', 'exon_index'}

    transcript_info: TranscriptInfo
    genomic_range: GenomicRange
    exon_index: int

    @property
    def transcript_id(self) -> str:
        return self.transcript_info.transcript_id

    @property
    def gene_id(self) -> str:
        return self.transcript_info.gene_id


@dataclass(init=False)
class CDSContextRepository:
    __slots__ = {'cds_ranges', '_target_ranges', '_target_cds_contexts'}

    cds_ranges: PyRanges
    _target_ranges: Optional[PyRanges]
    _target_cds_contexts: Optional[Dict[GenomicRange, Tuple[ExonInfo, GenomicRangePair]]]

    def __init__(self, cds_ranges: PyRanges) -> None:
        self.cds_ranges = cds_ranges
        self._target_ranges = None
        self._target_cds_contexts = None

    def register_target_ranges(self, target_ranges: PyRanges) -> None:
        self._target_ranges = target_ranges

    # TODO: investigate whether binding the chromosome and strand to each transcript
    # would speed up things (only one DataFrame to go through)
    def get_cds_by_index(self, transcript_id: str, exon_index: int) -> GenomicRange:
        exon_ranges: PyRanges = self.cds_ranges[
            (self.cds_ranges.transcript_id == transcript_id)
            & (self.cds_ranges.exon_index == exon_index)
        ]

        if not exon_ranges:
            raise RuntimeError(f"Exon index {exon_index} for transcript '{transcript_id}' not found!")

        rows: np.ndarray = exon_ranges.as_df()[[
            'Chromosome',
            'Strand',
            'Start',
            'End'
        ]].to_numpy()

        n = len(rows)

        if n > 1:
            raise RuntimeError(f"Exon index {exon_index} is not unique for transcript '{transcript_id}'!")

        row = rows[0]
        return GenomicRange(row[0], row[2] + 1, row[3], row[1])

    def get_cds_genomic_ranges(
        self,
        transcript_id: str,
        strand: str,
        exon_index: int,
        genomic_range: GenomicRange,
        delta_5p: int,
        delta_3p: int,
        cds_ext_5_length: int,
        cds_ext_3_length: int
    ) -> GenomicRangePair:

        def get_ext_5() -> Optional[GenomicRange]:

            def get_ext_5_from_previous_exon() -> GenomicRange:
                return self.get_cds_by_index(
                    transcript_id,
                    exon_index - 1
                ).get_from_3_prime(cds_ext_5_length)

            if cds_ext_5_length == 0:
                return None

            if exon_index == 0:
                raise ValueError("The first exon can't be out-of-frame!")

            if delta_5p == 1 and cds_ext_5_length > 1:
                raise ValueError("Unsupported partial exon: CDS extension would include both same and previous exon nucleotides!")

            return (
                get_ext_5_from_previous_exon() if delta_5p == 0 else
                genomic_range.get_before_5_prime(cds_ext_5_length)
            )

        def get_ext_3() -> Optional[GenomicRange]:

            def get_ext_3_from_next_exon() -> GenomicRange:
                return self.get_cds_by_index(
                    transcript_id,
                    exon_index + 1
                ).get_from_5_prime(cds_ext_3_length)

            if cds_ext_3_length == 0:
                return None

            if delta_3p == 1 and cds_ext_3_length > 1:
                raise ValueError("Unsupported partial exon: CDS extension would include both same and next exon nucleotides!")

            return (
                get_ext_3_from_next_exon() if delta_3p == 0 else
                genomic_range.get_past_3_prime(cds_ext_3_length)
            )

        validate_strand(strand)

        if cds_ext_5_length == 0 and cds_ext_3_length == 0:
            return None, None

        return (
            (get_ext_5(), get_ext_3()) if strand == '+' else
            (get_ext_3(), get_ext_5())
        )

    def compute_cds_contexts(self) -> None:

        def get_delta_start(df: pd.DataFrame, mask: pd.Series) -> pd.Series:
            return df.loc[mask, 'Start'] - df.loc[mask, 'Start_b']

        def get_delta_end(df: pd.DataFrame, mask: pd.Series) -> pd.Series:
            return df.loc[mask, 'End_b'] - df.loc[mask, 'End']

        if self._target_ranges is None:
            return

        exonic_ranges: pd.DataFrame = self._target_ranges.join(
            self.cds_ranges, strandedness='same').drop(['Strand_b', 'is_const']).as_df()

        # Skip if no CDS targets are found
        if exonic_ranges.shape[0] == 0:
            self._target_cds_contexts = {}
            return

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

        target_cds_extension_info: Dict[GenomicRange, Tuple[str, str, str, int, int, int, int, int]] = {
            pyrange_to_genomic_range(k, record): (
                record.gene_id,
                record.transcript_id,
                record.Strand,
                record.exon_index,
                record.delta_5p,
                record.delta_3p,
                record.frame,
                record.cds_ext_3_length
            )
            for k, ranges in PyRanges(df=exonic_ranges)
            for record in ranges.itertuples()
        }
        del exonic_ranges

        # TODO: optimise the retrieval of exons by transcript and index
        self._target_cds_contexts = {
            genomic_range: (
                ExonInfo(TranscriptInfo(gene_id, transcript_id), genomic_range, exon_index),
                self.get_cds_genomic_ranges(
                    transcript_id, strand, exon_index, genomic_range, delta_5p, delta_3p, cds_ext_5_length, cds_ext_3_length))
            for genomic_range, (
                gene_id,
                transcript_id,
                strand,
                exon_index,
                delta_5p,
                delta_3p,
                cds_ext_5_length,
                cds_ext_3_length
            ) in target_cds_extension_info.items()
        }

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


@dataclass
class UTRRepository:
    __slots__ = {'utr_ranges'}

    utr_ranges: PyRanges

    def get_transcript_infos(self, ref_ranges: PyRanges) -> Dict[GenomicRange, TranscriptInfo]:
        if not ref_ranges:
            return {}

        matches: PyRanges = ref_ranges.join(self.utr_ranges, strandedness='same')

        if not matches:
            return {}

        return {
            pyrange_to_genomic_range(k, record): TranscriptInfo(
                record.gene_id,
                record.transcript_id)
            for k, ranges in matches
            for record in ranges.itertuples()
        }


@dataclass
class AnnotationRepository:
    __slots__ = {'cds', 'utr'}

    cds: Optional[CDSContextRepository]
    utr: Optional[UTRRepository]

    @classmethod
    def from_ranges(
        cls,
        cds_ranges: Optional[PyRanges] = None,
        utr_ranges: Optional[PyRanges] = None
    ) -> AnnotationRepository:
        return cls(CDSContextRepository(cds_ranges), UTRRepository(utr_ranges))

    @classmethod
    def from_gff(cls, fp: str) -> AnnotationRepository:
        cds_ranges, utr_ranges = load_gff_cds(fp)
        return cls.from_ranges(cds_ranges, utr_ranges)
