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

from valiant.genomic_position_offsets import GenomicPositionOffsets, PosOffset, get_pos_offset
from valiant.uint_range import UIntRange
from valiant.var_stats import VarStats, get_alt_ref_delta


def get_del(pos, ref_len = 1):
    return VarStats(pos, ref_len, 0)


def get_ins(pos, alt_len = 1):
    return VarStats(pos, 0, alt_len)


pos_offsets = [PosOffset(10, -1), PosOffset(13, 0), PosOffset(15, 2)]


@pytest.mark.parametrize('pos,exp', [
    (9, 0),
    (10, -1),
    (11, -1),
    (14, 0),
    (15, 2),
    (20, 2)
])
def test_get_pos_offset(pos, exp):
    assert get_pos_offset(pos_offsets, pos) == exp


# XAY AAXAAY  AY AA
#  AXXAA AAXXXAXXAA
def test_genomic_position_offsets():
    var_stats = [
        get_del(10),
        get_ins(12),
        get_del(15),
        get_ins(18, 2),
        get_ins(20, 1)
    ]
    gpo = GenomicPositionOffsets.from_var_stats(var_stats, UIntRange(10, 22))
    print(gpo._ref_del_mask)
    print(gpo._pos_offsets)
    print(gpo._alt_ins_mask)

    print('\t'.join(['ref', 'alt', 'ref', 'alt']))
    for pos in range(gpo.ref_start, gpo.ref_range.end):
        alt = gpo.ref_to_alt_position(pos)
        ref = gpo.alt_to_ref_position(alt) if alt else None
        alt2 = gpo.ref_to_alt_position(ref) if ref else None
        if ref is not None:
            assert alt2 == alt
        print('\t'.join(str(x) if x is not None else '.' for x in [pos, alt, ref, alt2]))

    print('\t'.join(['alt', 'ref', 'alt']))
    for alt_pos in range(gpo.ref_start, gpo.ref_start + gpo.alt_length):
        ref = gpo.alt_to_ref_position(alt_pos)
        alt = gpo.ref_to_alt_position(ref) if ref is not None else None
        print('\t'.join(str(x) if x is not None else '.' for x in [alt_pos, ref, alt]))
        if alt is not None:
            assert alt == alt_pos

    assert len(gpo._ref_del_mask) == len(gpo.ref_range)
    assert len(gpo._alt_ins_mask) == len(gpo.ref_range) + get_alt_ref_delta(var_stats)

    # Position does not exist in the reference
    for pos in [11, 17, 18]:
        assert gpo.alt_to_ref_position(pos) is None

    # Position exists in the reference (unaltered)
    assert gpo.alt_to_ref_position(12) == 12
    assert gpo.ref_to_alt_position(12) == 12

    # Position exists in the reference (altered)
    assert gpo.alt_to_ref_position(19) == 18
    assert gpo.ref_to_alt_position(18) == 19


@pytest.mark.parametrize('variant,shift_range,del_range,ins_range', [
    (VarStats(15, 1, 1), None, None, None),  # SNV
    (VarStats(15, 5, 0), UIntRange(15, 19), UIntRange(15, 19), None),  # deletion
    (VarStats(15, 0, 5), UIntRange(15, 15), UIntRange(15, 15), UIntRange(15, 19)),  # insertion
    (VarStats(15, 5, 3), UIntRange(15, 19), UIntRange(15, 19), UIntRange(15, 19)),  # MNV (shift)
    (VarStats(15, 3, 5), UIntRange(15, 17), UIntRange(15, 17), UIntRange(15, 19)),  # MNV (shift)
    (VarStats(15, 5, 5), None, UIntRange(15, 19), UIntRange(15, 19))  # MNV (no shift)
])
def test_genomic_position_offsets_mask(variant, shift_range, del_range, ins_range):
    r = UIntRange(10, 20)
    gpo = GenomicPositionOffsets.from_var_stats([variant], r)

    for pos in r.positions:
        assert gpo._shift_mask[pos - gpo.ref_start] == (
            int(pos in shift_range) if shift_range is not None else 0)
        if variant.alt_ref_delta < 0:
            assert gpo._ref_del_mask[pos - gpo.ref_start] == int(pos in del_range)

    assert gpo.alt_range is not None
    if ins_range is not None and variant.alt_ref_delta > 0:
        for pos in gpo.alt_range.positions:
            assert gpo._alt_ins_mask[pos - gpo.ref_start] == int(pos in ins_range)
    else:
        for pos in gpo.alt_range.positions:
            assert gpo._alt_ins_mask[pos - gpo.ref_start] == 0
