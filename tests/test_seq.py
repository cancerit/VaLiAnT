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

from valiant.seq import Seq
from valiant.strings.dna_str import DnaStr
from valiant.uint_range import UIntRange


SEQ = 'ACGTACGT'


@pytest.mark.parametrize('start,end,exp', [
    (0, 1, 'AC'),
    (1, 3, 'CGT')
])
def test_seq_substr(start, end, exp):
    seq = Seq(0, DnaStr(SEQ))
    r = UIntRange(start, end)
    assert seq.substr(r) == exp
