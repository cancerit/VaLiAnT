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

import pytest
from contextlib import nullcontext
from valiant import utils


@pytest.mark.parametrize('strand,is_valid', [
    ('+', True),
    ('-', True),
    ('', False),
    (None, False),
    ('**', False)
])
def test_validate_strand(strand, is_valid):
    with pytest.raises(ValueError) if not is_valid else nullcontext():
        utils.validate_strand(strand)


@pytest.mark.parametrize('seq,rc', [
    ('AAA', 'TTT'),
    ('AAT', 'ATT'),
    ('ATT', 'AAT')
])
def test_reverse_complement(seq, rc):
    assert utils.reverse_complement(seq) == rc


@pytest.mark.parametrize('genomic_range,region', [
    (('X', 1, 100), 'X:1-100')
])
def test_get_region(genomic_range, region):
    assert utils.get_region(*genomic_range) == region


@pytest.mark.parametrize('seq,is_valid', [
    ('AAA', True),
    ('ACGT', True),
    ('acgt', False),
    (' ACGT', False),
    ('ACGT ', False),
    ('ACGTN', False),
    ('', False)
])
def test_is_dna(seq, is_valid):
    assert utils.is_dna(seq) == is_valid
