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

from valiant.genomic_position_offsets import GenomicPositionOffsets
from valiant.uint_range import UIntRange
from valiant.var_stats import VarStats


def get_del(pos, ref_len = 1):
    return VarStats(pos, ref_len, 0)


def get_ins(pos, alt_len = 1):
    return VarStats(pos, 0, alt_len)


def check_alt_to_ref_to_alt(gpo, ref):
    alt_1 = gpo.ref_to_alt_position(ref)
    print(f"{ref} -> {alt_1}")
    alt_2 = gpo.alt_to_ref_position(alt_1)
    print(f"{alt_1} -> {alt_2}")
    assert alt_1 == alt_2


def test_genomic_position_offsets():
    gpo = GenomicPositionOffsets.from_var_stats([
        get_del(10),
        get_ins(12),
        get_del(15),
        get_ins(18, 2),
        get_ins(20, 1)
    ], UIntRange(10, 22))
    assert gpo.alt_to_ref_position(11) is None
    assert gpo.alt_to_ref_position(12) == 12
    assert gpo.alt_to_ref_position(17) is None
    assert gpo.alt_to_ref_position(18) is None
    assert gpo.alt_to_ref_position(19) == 18
