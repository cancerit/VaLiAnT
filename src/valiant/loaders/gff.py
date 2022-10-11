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

import logging
from typing import Dict, Set, Tuple
import numpy as np
import pandas as pd
from pyranges import PyRanges, read_gtf
from ..utils import get_frame_complement, i8_0, i8_2, validate_strand

GFF_FIELDS = [
    'Chromosome',
    'Strand',
    'Start',
    'End',
    'Frame',
    'Feature',
    'gene_id',
    'transcript_id'
]

GFF_FEATURES = {
    'CDS',
    'stop_codon',
    'UTR'
}


def get_exon_indices(g: pd.core.groupby.generic.DataFrameGroupBy) -> pd.Series:

    # Validate strand information
    strand_counts: pd.Series = g.Strand.value_counts()

    n = strand_counts.shape[0]
    if n == 0:
        raise ValueError("No strand information available!")
    if n > 1:
        raise ValueError("Inconsistent strand information for transcript!")

    strand: str = strand_counts.index.get_level_values(1).to_numpy()[0]
    validate_strand(strand)

    # Assign indices based on strand
    is_plus_strand: bool = strand == '+'
    return g.cumcount(ascending=is_plus_strand).astype(np.int32)


def _get_last_cds_indices(cds_ranges: pd.DataFrame, strand: str) -> np.ndarray:
    return cds_ranges[cds_ranges.Strand == strand].groupby('transcript_id', sort=False).exon_index.idxmax()


def _get_frames(frame: pd.Series) -> pd.Series:
    # Convert frame
    try:
        frame = frame.astype(np.int8)
    except ValueError:
        raise ValueError("Invalid frame in GTF/GFF2 file: not an integer!")

    # Validate frame range
    if not frame.between(i8_0, i8_2, inclusive=True).all():
        raise ValueError("Invalid frame in GTF/GFF2 file: out of range!")

    # Transform frame to its complement
    return get_frame_complement(frame)


def load_gff_cds(fp: str) -> Tuple[PyRanges, PyRanges]:

    # Load necessary fields from GTF/GFF2 file
    ranges: pd.DataFrame = read_gtf(fp, as_df=True)[GFF_FIELDS].rename(columns={
        'Frame': 'frame'
    })

    logging.debug("GTF/GFF2 file: %d features found." % ranges.shape[0])

    # Drop unnecessary features
    ranges = ranges[ranges.Feature.isin(GFF_FEATURES)]

    logging.debug("GTF/GFF2 file: %d CDS features found." % ranges.shape[0])

    # Compress identifiers
    ranges.transcript_id = ranges.transcript_id.astype('category')
    ranges.gene_id = ranges.gene_id.astype('category')

    # Extract UTR features
    utr_mask: pd.Series = ranges.Feature == 'UTR'
    utr_ranges: PyRanges = PyRanges(df=ranges[utr_mask].drop('Feature', axis=1))
    ranges = ranges[~utr_mask]
    del utr_mask

    # Transform frames
    ranges.frame = _get_frames(ranges.frame)
    assert ranges.frame.dtype == 'int8'

    # Extract stop codon features
    stop_mask: pd.Series = ranges.Feature == 'stop_codon'
    stop_codons: Dict[str, int] = {
        r.transcript_id: r.End if r.Strand == '+' else r.Start
        for r in ranges[stop_mask].itertuples()
    }
    ranges = ranges[~stop_mask]
    del stop_mask

    # Sort CDS features by genomic coordinates
    ranges = ranges.drop('Feature', axis=1).sort_values(by=[
        'gene_id', 'transcript_id', 'Chromosome', 'Strand', 'Start', 'End'
    ], ignore_index=True)

    # Validate number of gene and transcript identifiers
    gene_n: int = ranges.gene_id.cat.categories.size
    transcript_n: int = ranges.transcript_id.cat.categories.size

    logging.debug("GTF/GFF2 file: %d genes found." % gene_n)
    logging.debug("GTF/GFF2 file: %d transcripts found." % transcript_n)

    if gene_n == 0 or transcript_n == 0:
        raise ValueError("No gene or transcript ID found in GTF/GFF file!")

    if transcript_n > gene_n:
        raise ValueError("Multiple transcripts per gene in GTF/GFF file are not supported!")

    # Check for missing identifiers
    if ranges.gene_id.isnull().values.any():
        raise ValueError("Missing gene ID in GTF/GFF2 file!")
    if ranges.transcript_id.isnull().values.any():
        raise ValueError("Missing transcript ID in GTF/GFF2 file!")

    # Check only one transcript per gene is present
    gene_transcript_counts: pd.DataFrame = ranges.groupby([
        'gene_id',
        'transcript_id'
    ], sort=False).size().reset_index(name='counts')
    gene_transcript_counts = gene_transcript_counts[gene_transcript_counts.counts > 0]
    if gene_transcript_counts.shape[0] > gene_n:
        raise ValueError("Multiple transcripts per gene in GTF/GFF file are not supported!")

    # Assign a sequential index to each CDS feature (5' to 3')
    ranges['exon_index'] = ranges.groupby(
        ['transcript_id'], sort=False).pipe(get_exon_indices)

    # Append the stop codons to the last CDS features
    strands: Set[str] = set(ranges.Strand.cat.categories.values)

    if '+' in strands:
        last_cds_plus: np.ndarray = _get_last_cds_indices(ranges, '+')
        ranges.loc[last_cds_plus, 'End'] += 3
        del last_cds_plus

    if '-' in strands:
        last_cds_minus: np.ndarray = _get_last_cds_indices(ranges, '-')
        ranges.loc[last_cds_minus, 'Start'] -= 3
        del last_cds_minus

    # TODO: validate with stop codon features
    del stop_codons

    # Convert to PyRanges
    cds_ranges: PyRanges = PyRanges(df=ranges)

    return cds_ranges, utr_ranges
