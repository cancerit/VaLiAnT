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
from valiant.loaders.gff import load_gff_cds
from .constants import GTF_SINGLE, GTF_MULTI
from .utils import get_data_file_path


@pytest.mark.parametrize('fp,valid', [
    (GTF_SINGLE, True),
    (GTF_MULTI, False)
])
def test_gff_load(fp, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():
        load_gff_cds(get_data_file_path(fp))
