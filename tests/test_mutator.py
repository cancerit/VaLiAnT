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

from valiant.mutators.snv_re import SnvReMutator
from valiant.codon_table import CodonTable
from valiant.codon_table_loader import load_codon_table_rows
from valiant.utils import get_default_codon_table_path
from valiant.seq import Seq
from valiant.strings.dna_str import DnaStr

codon_table = CodonTable.from_list(
    load_codon_table_rows(get_default_codon_table_path()))

seq = Seq(10, DnaStr('AAAAAA'))


def test_snv_re_mutator():
    mutator = SnvReMutator(codon_table)
    vars = mutator.get_variants(seq)
    for v in vars:
        assert v.ref_len == 3
        assert v.alt_len == 3
        assert (v.pos - seq.start) % 3 == 0
