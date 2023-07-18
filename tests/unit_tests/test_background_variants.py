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

from valiant.models.background_variants import _filter_variants_by_range, _get_alt_length, compute_genomic_offset, _compute_alt_offsets
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
