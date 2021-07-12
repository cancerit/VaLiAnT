########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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
from valiant.models.base import PositionRange, StrandedPositionRange


@pytest.mark.parametrize('start,end,is_valid', [
    (0, 12, False),
    (1, 12, True),
    (12, 1, False),
    (None, 1, False)
])
def test_position_range_init(start, end, is_valid):
    with pytest.raises((ValueError, TypeError)) if not is_valid else nullcontext():
        PositionRange(start, end)


@pytest.mark.parametrize('start,end,strand,is_valid', [
    (0, 12, '+', False),
    (1, 12, '+', True),
    (12, 1, '+', False),
    (None, 1, '+', False),
    (1, 12, '-', True),
    (1, 12, ':P', True),
    (1, 12, None, True)
])
def test_stranded_position_range_init(start, end, strand, is_valid):
    with pytest.raises((ValueError, TypeError)) if not is_valid else nullcontext():
        StrandedPositionRange(start, end, '+')


@pytest.mark.parametrize('start,end,strand', [
    (1, 10, '+')
])
def test_stranded_position_range_eq(start, end, strand):
    assert StrandedPositionRange(start, end, strand) == StrandedPositionRange(start, end, strand)
