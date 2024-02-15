########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

from valiant.strings.dna_str import DnaStr


s = DnaStr('AACC')

@pytest.mark.parametrize('n,exp', [
    (0, ''),
    (1, 'C'),
    (2, 'CC'),
    (4, 'AACC')
])
def test_tail(n, exp):
    assert s.tail(n) == exp


@pytest.mark.parametrize('n,exp', [
    (0, ''),
    (1, 'A'),
    (2, 'AA'),
    (4, 'AACC')
])
def test_head(n, exp):
    assert s.head(n) == exp


@pytest.mark.parametrize('n,exp', [
    ([0], 'A'),
    ([1, 2], 'AC'),
    ([3, 0, 2], 'CAC')
])
def test_get_at(n, exp):
    assert s.get_at(n) == exp
