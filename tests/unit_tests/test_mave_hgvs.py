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

import pytest
from valiant.enums import VariantType
from valiant.mave_hgvs import MAVEPrefix, get_mave_nt


@pytest.mark.parametrize('var_type,ref,alt,mave_nt_exp', [
    (VariantType.DELETION, 'AA', None, 'g.100_101del'),
    (VariantType.DELETION, 'A', None, 'g.100del'),
    (VariantType.SUBSTITUTION, 'A', 'T', 'g.100A>T'),
    (VariantType.SUBSTITUTION, 'A', 'TT', 'g.100delinsTT'),
    (VariantType.SUBSTITUTION, 'AA', 'TT', 'g.100_101delinsTT'),
    (VariantType.INSERTION, None, 'T', 'g.99_100insT'),
    (VariantType.INSERTION, None, 'TT', 'g.99_100insTT')
])
def test_get_mave_nt(var_type, ref, alt, mave_nt_exp):
    mave_nt = get_mave_nt(
        MAVEPrefix.LINEAR_GENOMIC, var_type.value, 100, ref, alt)

    assert mave_nt == mave_nt_exp
