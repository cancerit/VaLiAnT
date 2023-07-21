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

from valiant.models.base import GenomicPosition
from valiant.models.pam_protection import PamVariant
from valiant.models.variant import SubstitutionVariant
from valiant.models.new_pam import LAYER_BG, LAYER_PAM, _get_variant_groups


def test_alt_var_group_const():
    assert LAYER_BG == 0
    assert LAYER_PAM == 1


def test_get_variant_groups():
    bg_var = SubstitutionVariant(GenomicPosition('X', 100), 'A', 'T')
    pam_var = PamVariant.from_substitution(bg_var, 'sgrna-id-a')
    var_groups = _get_variant_groups([bg_var], [pam_var])
    assert var_groups[LAYER_BG].variants[0] == bg_var
    assert var_groups[LAYER_PAM].variants[0] == pam_var
