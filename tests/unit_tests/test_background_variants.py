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
from valiant.models.background_variants import GenomicPositionOffsets, _filter_variants_by_range, _get_alt_length, compute_genomic_offset, _compute_alt_offsets
from valiant.models.base import GenomicPosition
from valiant.models.variant import DeletionVariant, InsertionVariant


def get_pos(x):
    return GenomicPosition('X', x)


def test_compute_genomic_offset():
    assert compute_genomic_offset([
        DeletionVariant(get_pos(150), 'AA'),
        InsertionVariant(get_pos(100), 'AAA')
    ]) == 1


def test_compute_alt_offsets():
    variants = [
        DeletionVariant(get_pos(120), 'AA'),
        InsertionVariant(get_pos(110), 'AAA'),
        InsertionVariant(get_pos(130), 'AAA'),
        InsertionVariant(get_pos(999), 'AAA')
    ]
    ref_start = 100
    ref_length = 40
    variants_in_range = _filter_variants_by_range(
        ref_start, ref_length, variants, sort=True)
    alt_length = _get_alt_length(ref_length, variants_in_range)
    alt_offsets = _compute_alt_offsets(ref_start, alt_length, variants_in_range)
    assert alt_offsets.shape[0] == ref_length + compute_genomic_offset(variants_in_range)
    assert alt_offsets.tolist() == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]


def test_genomic_position_offsets_from_variants():
    variants = [
        DeletionVariant(get_pos(120), 'AA'),
        InsertionVariant(get_pos(110), 'AAA'),
        InsertionVariant(get_pos(130), 'AAA'),
        InsertionVariant(get_pos(999), 'AAA')
    ]
    ref_start = 100
    ref_length = 40
    gpo = GenomicPositionOffsets.from_variants(
        ref_start, ref_length, variants)

    assert gpo.ref_start == ref_start
    assert gpo.ref_length == ref_length
    assert gpo._alt_length == ref_length + 3 + 3 - 2
    assert len(gpo.variants_in_range) == 3
    assert len(gpo._pos_offsets) == len(gpo.variants_in_range)
    assert gpo._ins_offsets.shape[0] == gpo._alt_length


def test_genomic_position_offsets_alt_to_ref_position():
    variants = [
        InsertionVariant(get_pos(110), 'AAA'),
        DeletionVariant(get_pos(120), 'AA'),
        InsertionVariant(get_pos(130), 'AAA')
    ]
    ref_start = 100
    ref_length = 40
    gpo = GenomicPositionOffsets.from_variants(
        ref_start, ref_length, variants)

    positions = [
        variant.genomic_position.position
        for variant in variants
    ]

    with pytest.raises(ValueError):
        gpo.alt_to_ref_position(ref_start + gpo.alt_length)

    with pytest.raises(ValueError):
        gpo.alt_to_ref_position(ref_start - 1)

    for alt_pos, offset in [
        (ref_start, 0),
        (positions[0] + 3, 3),
        (positions[1] - 1, 3),
        (positions[2] + 3, 6),
    ]:
        assert gpo.alt_to_ref_position(alt_pos) == alt_pos - offset
