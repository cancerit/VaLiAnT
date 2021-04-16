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

import pytest
from valiant.enums import TargetonMutator
from valiant.models.refseq_ranges import ReferenceSequenceRanges, ReferenceSequenceRangeCollection


def test_reference_sequence_ranges_init():
    r2_mutators = {TargetonMutator.DEL1}
    mutators = set(), r2_mutators, set()
    chromosome = 'X'
    strand = '+'
    ref_start = 1
    ref_end = 1000
    r2_start = 100
    r2_end = 200
    c1_start = 1
    c1_end = 99
    c2_start = 221
    c2_end = 1000
    r3_start = 201
    r3_end = 220
    sgrna_ids = {'ID1', 'ID2'}

    # Initialise object
    rsr = ReferenceSequenceRanges(
        chromosome,
        strand,
        ref_start,
        ref_end,
        r2_start,
        r2_end,
        (0, 20),
        mutators,
        sgrna_ids)

    # Check reference range
    assert rsr.ref_range.chromosome == chromosome
    assert rsr.ref_range.strand == strand
    assert rsr.ref_range.start == ref_start
    assert rsr.ref_range.end == ref_end

    # Check sgRNA ID's
    assert rsr.sgrna_ids == sgrna_ids

    # Check constant regions
    assert len(rsr._const_regions) == 2
    c1, c2 = rsr._const_regions
    assert c1.chromosome == chromosome
    assert c1.strand == strand
    assert c1.start == c1_start
    assert c1.end == c1_end
    assert c2.chromosome == chromosome
    assert c2.strand == strand
    assert c2.start == c2_start
    assert c2.end == c2_end

    # Check target regions
    assert len(rsr._target_regions) == 3
    r1, r2, r3 = rsr._target_regions
    assert r1 is None
    assert r2.mutators == r2_mutators
    assert r2.genomic_range.chromosome == chromosome
    assert r2.genomic_range.strand == strand
    assert r2.genomic_range.start == r2_start
    assert r2.genomic_range.end == r2_end
    assert not r3.mutators
    assert r3.genomic_range.chromosome == chromosome
    assert r3.genomic_range.strand == strand
    assert r3.genomic_range.start == r3_start
    assert r3.genomic_range.end == r3_end


@pytest.mark.parametrize('s,exp_mutators', [
    ('(1del), (1del,snv), ()', [
        {TargetonMutator.DEL1},
        {TargetonMutator.DEL1, TargetonMutator.SNV},
        set()
    ])
])
def test_reference_sequence_ranges_parse_mutators(s, exp_mutators):
    mutators = ReferenceSequenceRanges.parse_mutator_tuples(s)
    assert mutators == exp_mutators


def test_reference_sequence_range_collection_init():
    r2_mutators = {TargetonMutator.DEL1}
    mutators = set(), r2_mutators, set()
    chromosome = 'X'
    strand = '+'
    ref_start = 1
    ref_end = 1000
    r2_start = 100
    r2_end = 200
    sgrna_ids = {'ID1', 'ID2'}

    rsr = ReferenceSequenceRanges(
        chromosome,
        strand,
        ref_start,
        ref_end,
        r2_start,
        r2_end,
        (0, 20),
        mutators,
        sgrna_ids)

    rsrc = ReferenceSequenceRangeCollection([rsr])

    # Check reference sequence ranges
    assert len(rsrc._rsrs) == 1
    assert rsrc._rsrs[0] == rsr

    # Check sgRNA ID's
    assert rsrc.sgrna_ids == sgrna_ids

    # Check range tables
    assert len(rsrc._ref_ranges) == 1
    assert len(rsrc._region_ranges) == 4
