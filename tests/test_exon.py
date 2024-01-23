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

from valiant.exon import Exon
from valiant.strings.strand import Strand
from valiant.uint_range import UIntRange


plus = Strand('+')
minus = Strand('-')


@pytest.mark.parametrize('strand,codon_index,exp', [
    (plus, 0, UIntRange(10, 12)),
    (minus, 0, UIntRange(14, 16)),
    (plus, 1, UIntRange(13, 15)),
    (minus, 1, UIntRange(11, 13)),
    # Partial exon
    (plus, 2, UIntRange(16, 18))
])
def test_exon_get_codon(strand, codon_index, exp):
    exon = Exon(10, 16, 0, 0)
    assert exon.get_codon(strand, codon_index) == exp
