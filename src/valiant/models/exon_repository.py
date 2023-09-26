from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from pyranges import PyRanges

from .base import GenomicRange, GenomicRangePair
from .exon_ext_info import ExonExtInfo
from ..utils import validate_strand


@dataclass(frozen=True)
class ExonRepository:
    __slots__ = ['cds_ranges']

    cds_ranges: PyRanges

    # TODO: investigate whether binding the chromosome and strand to each transcript
    #  would speed up things (only one DataFrame to go through)
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

    def _get_ext_5(self, exon: ExonExtInfo) -> Optional[GenomicRange]:

        def get_ext_5_from_previous_exon() -> GenomicRange:
            return self.get_cds_by_index(
                exon.exon_info.transcript_id,
                exon.exon_info.exon_index - 1
            ).get_from_3_prime(exon.cds_ext_5_length)

        if exon.cds_ext_5_length == 0:
            return None

        if exon.exon_info.exon_index == 0:
            raise ValueError("The first exon can't be out-of-frame!")

        if exon.delta_5p == 1 and exon.cds_ext_5_length > 1:
            raise ValueError(
                "Unsupported partial exon: CDS extension would include both same and previous exon nucleotides!")

        return (
            get_ext_5_from_previous_exon() if exon.delta_5p == 0 else
            exon.exon_info.genomic_range.get_before_5_prime(exon.cds_ext_5_length)
        )

    def _get_ext_3(self, exon: ExonExtInfo) -> Optional[GenomicRange]:

        def get_ext_3_from_next_exon() -> GenomicRange:
            return self.get_cds_by_index(
                exon.exon_info.transcript_id,
                exon.exon_info.exon_index + 1
            ).get_from_5_prime(exon.cds_ext_3_length)

        if exon.cds_ext_3_length == 0:
            return None

        if exon.delta_3p == 1 and exon.cds_ext_3_length > 1:
            raise ValueError(
                "Unsupported partial exon: CDS extension would include both same and next exon nucleotides!")

        return (
            get_ext_3_from_next_exon() if exon.delta_3p == 0 else
            exon.exon_info.genomic_range.get_past_3_prime(exon.cds_ext_3_length)
        )

    def get_exon_ext_genomic_ranges(
        self,
        exon: ExonExtInfo
    ) -> GenomicRangePair:
        """
        Retrieve the genomic ranges of the non-adjacent nucleotides of the codons
        crossing the exon boundaries, if any
        """

        validate_strand(exon.strand)

        if exon.cds_ext_5_length == 0 and exon.cds_ext_3_length == 0:
            return None, None

        ext_5 = self._get_ext_5(exon)
        ext_3 = self._get_ext_3(exon)

        return (
            (ext_5, ext_3) if exon.strand == '+' else
            (ext_3, ext_5)
        )
