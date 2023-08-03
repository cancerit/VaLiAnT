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
from valiant.models.base import GenomicPosition
from valiant.models.new_pam import PamBgAltSeqBuilder
from valiant.models.pam_protection import PamVariant
from valiant.models.variant import InsertionVariant, DeletionVariant, SubstitutionVariant
from .constants import DUMMY_PAM_PROTECTION_NT
from .utils import get_ref_seq


SEQ = 'AAAACGTACGTACGT'

n = len(SEQ)
gp = GenomicPosition('X', 1)


@pytest.mark.parametrize('variant,ref,pam', [
    (InsertionVariant(gp, 'TTT'), 'TTT' + SEQ, 'TTT' + n * DUMMY_PAM_PROTECTION_NT),
    (DeletionVariant(gp, 'AAA'), SEQ[3:], (n - 3) * DUMMY_PAM_PROTECTION_NT),
    (SubstitutionVariant(gp, 'A', 'G'), 'G' + SEQ[1:], 'G' + (n - 1) * DUMMY_PAM_PROTECTION_NT)
])
@pytest.mark.parametrize('pam_protection', [True, False])
def test_apply_variant(variant, ref, pam, pam_protection):
    seq = SEQ
    start = 1
    ref_seq = get_ref_seq(seq, pos=start)
    pam_ref_seq = PamBgAltSeqBuilder.from_ref_seq(ref_seq, [], [
        PamVariant(GenomicPosition('X', start + i), seq[i], DUMMY_PAM_PROTECTION_NT, 'sgrna-id-a')
        for i in range(len(seq))
    ] if pam_protection else [])
    mseq = pam_ref_seq.mutate(variant)
    assert mseq == (pam if pam_protection else ref)
