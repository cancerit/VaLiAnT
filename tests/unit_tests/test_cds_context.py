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

import pandas as pd
from pyranges import PyRanges
import pytest

from valiant.models.base import GenomicRange, TranscriptInfo
from valiant.models.cds_context_repository import CDSContextRepository
from valiant.models.exon_info import ExonInfo
from valiant.models.exon_ext_info import ExonExtInfo
from valiant.models.exon_repository import ExonRepository

PYRANGES_FIELDS = ['Chromosome', 'Strand', 'Start', 'End']
GID = 'GENE_ID_001'
TID = 'TRANSCRIPT_ID_001'

RANGE = ('X', '+', 100, 120)
GR = GenomicRange('X', 100, 120, '+')

RANGES = pd.DataFrame.from_records([
    RANGE
], columns=PYRANGES_FIELDS)

CDS_RANGES = [
    ('X', '+', 100, 120, GID, TID, 0, 0),
    ('X', '+', 200, 207, GID, TID, 2, 1),
    ('X', '+', 300, 304, GID, TID, 1, 2)
]

CDS_RANGES_DF = pd.DataFrame.from_records([
    (chromosome, strand, start - 1, end, gene_id, transcript_id, frame, exon_index)
    for chromosome, strand, start, end, gene_id, transcript_id, frame, exon_index in CDS_RANGES
], columns=PYRANGES_FIELDS + ['gene_id', 'transcript_id', 'frame', 'exon_index'])


def test_exon_repository_init():
    ranges = PyRanges(df=RANGES)
    ccr = ExonRepository(ranges)

    assert ccr.cds_ranges == ranges


def test_exon_repository_get_cds_by_index():
    chromosome, strand, start, end, _, transcript_id, _, exon_index = CDS_RANGES[0]
    ccr = ExonRepository(PyRanges(df=CDS_RANGES_DF))
    gr = ccr.get_cds_by_index(transcript_id, exon_index)

    assert gr.chromosome == chromosome
    assert gr.strand == strand
    assert gr.start == start
    assert gr.end == end

    with pytest.raises(Exception):
        ccr.get_cds_by_index(transcript_id, 999)


@pytest.mark.parametrize('strand,len5p,len3p,exp_cds_pre,exp_cds_suf', [
    ('+', 0, 0, None, None)
])
def test_exon_repository_get_cds_genomic_ranges(strand, len5p, len3p, exp_cds_pre, exp_cds_suf):
    transcript_id = TID
    exon_index = 1
    exons = ExonRepository(PyRanges(df=CDS_RANGES_DF))
    cds_pre, cds_suf = exons.get_exon_ext_genomic_ranges(ExonExtInfo(
        TranscriptInfo("GENE_ID", transcript_id),
        GR,
        exon_index,
        len5p, len3p))

    assert cds_pre == exp_cds_pre
    assert cds_suf == exp_cds_suf


@pytest.mark.parametrize('start,end,exp_ext_5,exp_ext_3', [
    (100, 117, None, None),
    (100, 118, None, (119, 120)),
    (100, 119, None, (120, 120))
])
def test_cds_context_repository_compute_cds_contexts(start, end, exp_ext_5, exp_ext_3):
    cds_ranges = PyRanges(df=CDS_RANGES_DF)

    chromosome = 'X'
    strand = '+'
    gr = GenomicRange(chromosome, start, end, strand)
    target_ranges = PyRanges(df=pd.DataFrame.from_records([
        gr.as_pyrange()
    ], columns=PYRANGES_FIELDS))
    target_ranges.is_const = False

    # Initialise exon repository
    exons = ExonRepository(cds_ranges)

    # Compute CDS contexts
    ccr = CDSContextRepository.from_exons(target_ranges, exons)

    # Check CDS contexts
    assert len(ccr._target_cds_contexts) == 1
    exon_info, (ext_5, ext_3) = ccr._target_cds_contexts[gr]

    # Check exon information
    assert isinstance(exon_info, ExonInfo)
    assert exon_info.gene_id == GID
    assert exon_info.transcript_id == TID

    # Check CDS extension
    if exp_ext_5 is not None:
        assert ext_5 == GenomicRange(chromosome, *exp_ext_5, strand)
    else:
        assert ext_5 is None

    if exp_ext_3 is not None:
        assert ext_3 == GenomicRange(chromosome, *exp_ext_3, strand)
    else:
        assert ext_3 is None

    # Check information retrieval
    assert ccr.get_cds_extensions(gr) == (ext_5, ext_3)
    assert ccr.get_exon_info(gr) == exon_info
    assert ccr.get_transcript_info(gr) == exon_info.transcript_info
