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
import os
import pytest
from valiant.loaders.codon_table import _parse_codon_table_row, load_codon_table
from valiant.models.codon_table import CodonTable
from .constants import CODON_TABLE_FP
from .utils import get_data_file_path


@pytest.mark.parametrize('row,valid', [
    (('ACC', 'X', '0.23', 'RANK2'), True),
    (('TTT', 'STOP', '0.23', 'RANK2'), True),
    (('ACC', 'X', '0.23', 'RANKZ'), False),
    (('ACC', 'X', '2.23', 'RANK2'), False),
    (('ACC', 'XYZ', '0.23', 'RANK2'), False),
    (('ZZZ', 'X', '0.23', 'RANK2'), False)
])
def test_parse_codon_table_row(row, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():
        _parse_codon_table_row(row)


@pytest.mark.parametrize('fp,exists', [
    (CODON_TABLE_FP, True),
    ('missing.csv', False)
])
def test_load_codon_table(fp, exists):
    with pytest.raises(FileNotFoundError) if not exists else nullcontext():
        load_codon_table(get_data_file_path(fp))


@pytest.mark.parametrize('fp', [CODON_TABLE_FP])
def test_codon_table_init(fp):
    ct_rows = load_codon_table(get_data_file_path(fp))
    CodonTable(ct_rows)


@pytest.mark.parametrize('fp', [CODON_TABLE_FP])
def test_codon_table_load(fp):
    CodonTable.load(get_data_file_path(fp))
