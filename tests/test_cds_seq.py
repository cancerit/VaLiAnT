########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from valiant.cds_seq import CdsSeq
from valiant.strings.codon import Codon
from valiant.strings.dna_str import DnaStr


@pytest.mark.parametrize('seq,prefix,suffix,valid', [
    ('ACG', None, None, True),
    ('CG', 'A', None, True),
    ('C', 'A', 'G', True),
    ('A', None, 'CG', True),
    ('ACG', 'A', 'CG', False),
    ('AC', 'CG', 'GC', False),
    ('ACGTAC', None, None, False)
])
def test_cds_seq_as_codon(seq, prefix, suffix, valid):
    seq = CdsSeq(
        10,
        DnaStr(seq),
        cds_prefix=prefix or DnaStr.empty(),
        cds_suffix=suffix or DnaStr.empty())

    with pytest.raises(ValueError) if not valid else nullcontext():
        assert seq.as_codon() == Codon('ACG')
