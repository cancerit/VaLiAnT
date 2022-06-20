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

import numpy as np
import pytest
from valiant.enums import VariantType
from valiant.models.base import GenomicRange
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.targeton import Targeton, CDSTargeton
from .utils import get_no_op_pam_protected_sequence

del_var_type = np.int8(VariantType.DELETION.value)

del_offset_method = {
    0: 'get_2del0_mutations',
    1: 'get_2del1_mutations'
}


@pytest.mark.parametrize('offset,seq,exp_pos,exp_ref,exp_mseq', [
    (0, 'AAATTT', [0, 2, 4], ['AA', 'AT', 'TT'], ['ATTT', 'AATT', 'AAAT']),
    (1, 'AAATTT', [1, 3], ['AA', 'TT'], ['ATTT', 'AAAT']),
    (0, 'AATTT', [0, 2], ['AA', 'TT'], ['TTT', 'AAT']),
    (1, 'AATTT', [1, 3], ['AT', 'TT'], ['ATT', 'AAT'])
])
@pytest.mark.parametrize('cds', [True, False])
def test_get_2del_mutations(offset, seq, exp_pos, exp_ref, exp_mseq, cds):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, '+')
    pam_seq = get_no_op_pam_protected_sequence(seq, gr)
    t = (
        CDSTargeton.from_pam_seq(pam_seq, 'AA', 'A') if cds else
        Targeton.from_pam_seq(pam_seq)
    )

    # Generate in-frame deletions
    mc = getattr(t, del_offset_method[offset])()

    # Check metadata table
    assert np.array_equal(mc.df.mut_position, np.array(exp_pos))
    assert np.array_equal(mc.df.ref.astype('string').to_numpy(), np.array(exp_ref))
    assert np.array_equal(mc.df.mseq.astype('string').to_numpy(), np.array(exp_mseq))
    assert mc.df.new.isna().all()
    assert mc.df.var_type.unique()[0] == del_var_type


@pytest.mark.parametrize('seq,pre,suf,strand,exp_pos,exp_ref,exp_mseq', [
    ('AAATTT', '', '', '+', [0, 3], ['AAA', 'TTT'], ['TTT', 'AAA']),
    ('GAAATTT', 'GG', '', '+', [1, 4], ['AAA', 'TTT'], ['GTTT', 'GAAA']),
    ('AAATTTGG', '', 'G', '+', [0, 3], ['AAA', 'TTT'], ['TTTGG', 'AAAGG']),
    ('GAAATTT', 'GG', '', '-', [1, 4], ['AAA', 'TTT'], ['GTTT', 'GAAA']),
    ('AAATTTGG', '', 'G', '-', [0, 3], ['AAA', 'TTT'], ['TTTGG', 'AAAGG'])
])
def test_get_inframe_mutations(seq, pre, suf, strand, exp_pos, exp_ref, exp_mseq):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, strand)
    pam_seq = get_no_op_pam_protected_sequence(seq, gr)
    t = CDSTargeton.from_pam_seq(pam_seq, pre, suf)

    # Generate in-frame deletions
    mc = t.get_inframe_mutations()

    # Check metadata table
    assert np.array_equal(mc.df.mut_position, np.array(exp_pos))
    assert np.array_equal(mc.df.ref.astype('string').to_numpy(), np.array(exp_ref))
    assert np.array_equal(mc.df.mseq.astype('string').to_numpy(), np.array(exp_mseq))
    assert mc.df.new.isna().all()
    assert mc.df.var_type.unique()[0] == del_var_type
