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
