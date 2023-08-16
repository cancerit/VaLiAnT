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

from contextlib import nullcontext

import pytest
from valiant.models.uint_range import UIntRange


def test_uint_range_sort():
    a = UIntRange(10, 12)
    b = UIntRange(1, 8)
    assert sorted([a, b]) == [b, a]


def test_uint_range_diff():
    a = UIntRange(10, 100)

    x = UIntRange(10, 19)
    y = UIntRange(91, 100)
    z = UIntRange(51, 59)

    assert a.diff([]) == [a]

    assert a.diff([x, y, z]) == [
        UIntRange(20, 50),
        UIntRange(60, 90)
    ]


@pytest.mark.parametrize('skip', [True, False])
def test_uint_range_diff_skip(skip):
    a = UIntRange(100, 200)
    b = UIntRange(10, 20)

    with nullcontext() if skip else pytest.raises(ValueError):
        assert a.diff([b], skip_out_of_range=skip) == [a]
