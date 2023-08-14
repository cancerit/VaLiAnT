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

    assert a.diff([x, y, z]) == [
        UIntRange(20, 50),
        UIntRange(60, 90)
    ]
