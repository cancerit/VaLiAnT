########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from valiant.models.base import GenomicPosition


@pytest.mark.parametrize('chromosome,position,valid', [
    ('X', 100, True),
    ('', 100, False),
    ('X', 0, False),
    ('X', -100, False)
])
def test_genomic_position_init(chromosome, position, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():
        gp = GenomicPosition(chromosome, position)
        assert gp.chromosome == chromosome
        assert gp.position == position


def test_genomic_position_str():
    assert str(GenomicPosition('X', 10)) == 'X:10'
