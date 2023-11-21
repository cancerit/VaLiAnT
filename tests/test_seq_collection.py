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

from valiant.seq import Seq
from valiant.seq_collection import SeqCollection
from valiant.strings.dna_str import DnaStr
from valiant.uint_range import UIntRange


sc = SeqCollection([
    Seq(10, DnaStr('AC')),
    Seq(20, DnaStr('GGGGGG')),
    Seq(30, DnaStr('CA'))
])


def test_seq_collection_get_before():
    r = UIntRange(22, 24)
    # Test within the sequence boundary
    assert sc.get_before(1, r, 2) == DnaStr('GG')
    # Test across the sequence boundary
    assert sc.get_before(1, r, 3) == DnaStr('CGG')
    r = UIntRange(20, 24)
    # Test before the sequence boundary
    assert sc.get_before(1, r, 2) == DnaStr('AC')


def test_seq_collection_get_after():
    r = UIntRange(20, 23)
    # Test within the sequence boundary
    assert sc.get_after(1, r, 2) == DnaStr('GG')
    # Test across the sequence boundary
    assert sc.get_after(1, r, 3) == DnaStr('GGC')
    r = UIntRange(24, 25)
    # Test after the sequence boundary
    assert sc.get_after(1, r, 2) == DnaStr('CA')
