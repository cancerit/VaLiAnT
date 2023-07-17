########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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
from contextlib import nullcontext
from valiant.models.base import GenomicPosition, GenomicRange
from valiant.models.sequences import ReferenceSequence
from valiant.models.pam_protection import PamVariant, compute_pam_protected_sequence


@pytest.mark.parametrize('seq,pos,ref,alt,ppseq,valid', [
    ('AACAA', 102, 'C', 'T', 'AATAA', True),
    ('AACAA', 102, 'G', 'T', 'AATAA', False),  # unexpected reference nucleotide
    ('AACAA', 1000, 'C', 'T', 'AATAA', False)  # variant out of range
])
def test_compute_pam_protected_sequence(seq, pos, ref, alt, ppseq, valid):
    chromosome = 'X'
    start = 100
    end = start + len(seq) - 1

    variant = PamVariant(GenomicPosition(chromosome, pos), ref, alt, 'sgrna-1')

    gr = GenomicRange(chromosome, start, end, '+')
    ref_seq = ReferenceSequence(seq, gr)
    with pytest.raises(Exception) if not valid else nullcontext():
        pam_ref_seq = compute_pam_protected_sequence(ref_seq, {variant}, set())

    if valid:
        assert pam_ref_seq.sequence == ref_seq.sequence
        assert pam_ref_seq.pam_protected_sequence == ppseq
