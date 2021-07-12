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

import pytest
from valiant.models.codon_table import STOP_CODE
from valiant.snv import build_aa_sub_table, get_fast_aa_subs
from .utils import load_codon_table, get_aux_tables

codon_table = load_codon_table()


@pytest.mark.parametrize('strand', ['+', '-'])
def test_build_aa_sub_table(strand):
    t = build_aa_sub_table(codon_table, strand)
    assert t.shape[0] == 1219
    t.ref_aa.cat.remove_unused_categories(inplace=True)
    t.alt_aa.cat.remove_unused_categories(inplace=True)
    assert STOP_CODE in t.ref_aa.cat.categories
    assert STOP_CODE not in t.alt_aa.cat.categories


aux = get_aux_tables()
all_aa_table = aux.all_aa_table


@pytest.mark.parametrize('seq', ['AAACCCTGA'])
@pytest.mark.parametrize('frame', [0, 1, 2])
@pytest.mark.parametrize('prefix', [0, 1, 2])
@pytest.mark.parametrize('suffix', [0, 1, 2])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_get_fast_aa_subs(seq, frame, prefix, suffix, strand):
    seq_ = ('T' * prefix) + seq + ('T' * suffix)
    df = get_fast_aa_subs(
        all_aa_table._strand_table[strand], frame, seq_)
    assert not df.ref.isna().any()
