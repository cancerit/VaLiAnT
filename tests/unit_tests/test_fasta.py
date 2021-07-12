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
from valiant.errors import SequenceNotFound
from valiant.loaders.fasta import load_from_multi_fasta
from .constants import MULTI_FASTA_FP
from .utils import get_data_file_path


SEQ_ID_1 = 'seq001'
SEQ_ID_2 = 'seq002'

multi_fasta_fp = get_data_file_path(MULTI_FASTA_FP)


@pytest.mark.parametrize('seq_ids,all_exist', [
    ([SEQ_ID_1, SEQ_ID_2], True),
    ([SEQ_ID_1, 'missing-id'], False)
])
def test_load_from_multi_fasta(seq_ids, all_exist):
    with pytest.raises(SequenceNotFound) if not all_exist else nullcontext():
        id_seqs = load_from_multi_fasta(multi_fasta_fp, seq_ids)

    if all_exist:
        assert set(seq_ids) == set(id_seqs.keys())
