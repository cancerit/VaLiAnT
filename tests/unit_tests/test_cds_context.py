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

import pandas as pd
from pyranges import PyRanges
import pytest
from valiant.models.base import GenomicRange
from valiant.models.exon import ExonInfo, CDSContextRepository

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


def test_cds_context_repository_init():
    ranges = PyRanges(df=RANGES)
    ccr = CDSContextRepository(ranges)

    assert ccr.cds_ranges == ranges


def test_cds_context_repository_register_target_ranges():
    ranges = PyRanges(df=RANGES)
    ccr = CDSContextRepository(ranges)
    ccr.register_target_ranges(ranges)

    assert ccr._target_ranges == ranges


def test_cds_context_repository_get_cds_by_index():
    chromosome, strand, start, end, _, transcript_id, _, exon_index = CDS_RANGES[0]
    ccr = CDSContextRepository(PyRanges(df=CDS_RANGES_DF))
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
def test_cds_context_repository_get_cds_genomic_ranges(strand, len5p, len3p, exp_cds_pre, exp_cds_suf):
    transcript_id = TID
    exon_index = 1
    ccr = CDSContextRepository(PyRanges(df=CDS_RANGES_DF))
    cds_pre, cds_suf = ccr.get_cds_genomic_ranges(transcript_id, strand, exon_index, GR, 0, 0, len5p, len3p)

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

    # Initialise repository
    ccr = CDSContextRepository(cds_ranges)
    ccr.register_target_ranges(target_ranges)

    # Compute CDS contexts
    ccr.compute_cds_contexts()

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
