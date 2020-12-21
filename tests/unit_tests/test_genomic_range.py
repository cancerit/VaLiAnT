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

from contextlib import nullcontext
import pytest
from valiant.models.base import GenomicRange, GenomicPosition


@pytest.mark.parametrize('chromosome,start,end,strand,valid', [
    ('X', 1, 10, '+', True),
    ('X', 1, 10, '-', True),
    ('', 1, 10, '+', False),
    ('X', -10, 10, '+', False),
    ('X', 1, 10, '*', False)
])
def test_genomic_range_init(chromosome, start, end, strand, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():
        GenomicRange(chromosome, start, end, strand)


def test_genomic_range_region():
    assert GenomicRange('X', 1, 10, '+').region == 'X:1-10'


def test_genomic_range_length():
    assert len(GenomicRange('X', 5, 10, '+')) == 6


@pytest.mark.parametrize('a_range, b_range, b_contains_a', [
    ((100, 1000), (90, 1010), True),
    ((100, 1000), (90, 150), False)
])
def test_genomic_range_contains(a_range, b_range, b_contains_a):
    chromosome = 'X'
    strand = '+'
    a = GenomicRange(chromosome, *a_range, strand)
    b = GenomicRange(chromosome, *b_range, strand)
    assert a in a
    assert (a in b if b_contains_a else a not in b)
    assert b not in a


@pytest.mark.parametrize('n,exp_start,exp_end', [
    (3, 100, 102)
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_genomic_range_get_head(n, exp_start, exp_end, strand):
    chromosome = 'X'
    gr = GenomicRange(chromosome, 100, 200, strand)
    assert gr.get_head(n) == GenomicRange(
        chromosome, exp_start, exp_end, strand)


@pytest.mark.parametrize('n,exp_start,exp_end', [
    (3, 198, 200)
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_genomic_range_get_tail(n, exp_start, exp_end, strand):
    chromosome = 'X'
    gr = GenomicRange(chromosome, 100, 200, strand)
    assert gr.get_tail(n) == GenomicRange(
        chromosome, exp_start, exp_end, strand)


@pytest.mark.parametrize('strand,n,exp_start,exp_end', [
    ('+', 3, 100, 102),
    ('-', 3, 198, 200)
])
def test_genomic_range_get_from_5_prime(strand, n, exp_start, exp_end):
    chromosome = 'X'
    gr = GenomicRange(chromosome, 100, 200, strand)
    assert gr.get_from_5_prime(n) == GenomicRange(
        chromosome, exp_start, exp_end, strand)


@pytest.mark.parametrize('strand,n,exp_start,exp_end', [
    ('-', 3, 100, 102),
    ('+', 3, 198, 200)
])
def test_genomic_range_get_from_3_prime(strand, n, exp_start, exp_end):
    chromosome = 'X'
    gr = GenomicRange(chromosome, 100, 200, strand)
    assert gr.get_from_3_prime(n) == GenomicRange(
        chromosome, exp_start, exp_end, strand)


@pytest.mark.parametrize('parent,child,exp_range,valid', [
    (GenomicRange('X', 100, 200, '+'), GenomicRange('X', 150, 160, '+'), (50, 61), True),
    (GenomicRange('X', 100, 200, '+'), GenomicRange('X', 20, 160, '+'), None, False),
    (GenomicRange('X', 100, 200, '+'), GenomicRange('X', 150, 160, '-'), None, False)

])
def test_genomic_range_get_relative_subrange(parent, child, exp_range, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():
        assert exp_range == parent.get_relative_subrange(child)


@pytest.mark.parametrize('pos,exp', [
    (GenomicPosition('X', 100), True),
    (GenomicPosition('X', 200), True),
    (GenomicPosition('X', 150), True),
    (GenomicPosition('X', 300), False),
    (GenomicPosition('Y', 100), False)
])
def test_genomic_range_contains_position(pos,exp):
    gr = GenomicRange('X', 100, 200, '+')
    assert gr.contains_position(pos) == exp


@pytest.mark.parametrize('chromosome,start,end', [
    ('X', 100, 200)
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_genomic_range_as_unstranded(chromosome, start, end, strand):
    assert GenomicRange(chromosome, start, end, strand).as_unstranded() == (chromosome, start, end)


@pytest.mark.parametrize('chromosome,start,end,exp_start', [
    ('X', 100, 200, 99)
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_genomic_range_as_pyrange(chromosome, start, end, exp_start, strand):
    assert GenomicRange(chromosome, start, end, strand).as_pyrange() == (chromosome, strand, exp_start, end)


@pytest.mark.parametrize('chromosome,start,end,exp_start', [
    ('X', 100, 200, 99)
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_genomic_range_as_unstranded_pyrange(chromosome, start, end, exp_start, strand):
    assert GenomicRange(chromosome, start, end, strand).as_unstranded_pyrange() == (chromosome, exp_start, end)
