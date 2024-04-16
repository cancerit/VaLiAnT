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

import pytest

from valiant.seq import Seq
from valiant.seq_converter import apply_variants
from valiant.strings.dna_str import DnaStr
from valiant.variant import Variant


seq = Seq(10, DnaStr('A' * 10))


def get_del(pos, length):
	return Variant(pos, DnaStr('A' * length), DnaStr.empty())


@pytest.mark.parametrize('variants', [
	[get_del(10, 10)],
	[get_del(10, 5), get_del(15, 3)]
])
def test_apply_variants(variants):
	delta = sum(x.alt_ref_delta for x in variants)
	alt_length = len(seq) + delta
	alt_seq = apply_variants(seq, alt_length, variants)
	assert len(alt_seq) == alt_length
