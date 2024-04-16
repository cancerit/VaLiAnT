########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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

from valiant.strings.dna_str import DnaStr
from valiant.seq import Seq
from valiant.uint_range import UIntRange, UIntRangeSortedList


sequences = [
    Seq(10, DnaStr('AC')),
    Seq(20, DnaStr('GGGGGG')),
    Seq(30, DnaStr('CA'))
]


sc = UIntRangeSortedList([
    seq.get_range()
    for seq in sequences
])


def test_uint_range_sorted_list_get_before():
    r = UIntRange(22, 24)
    # Test within the sequence boundary
    assert sc.get_before(1, r, 2) == [20, 21]
    # Test across the sequence boundary
    assert sc.get_before(1, r, 3) == [11, 20, 21]
    r = UIntRange(20, 24)
    # Test before the sequence boundary
    assert sc.get_before(1, r, 2) == [10, 11]


def test_uint_range_sorted_list_get_after():
    r = UIntRange(20, 23)
    # Test within the sequence boundary
    assert sc.get_after(1, r, 2) == [24, 25]
    # Test across the sequence boundary
    assert sc.get_after(1, r, 3) == [24, 25, 30]
    r = UIntRange(24, 25)
    # Test after the sequence boundary
    assert sc.get_after(1, r, 2) == [30, 31]
