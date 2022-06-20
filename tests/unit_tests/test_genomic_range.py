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
def test_genomic_range_contains_position(pos, exp):
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
