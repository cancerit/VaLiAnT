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
from valiant.models.new_pam import PamBgAltSeqBuilder
from valiant.models.sequences import ReferenceSequence
from valiant.models.pam_protected_reference_sequence import PamVariant
from .utils import get_ref_seq


@pytest.mark.parametrize('seq,pos,ref,alt,ppseq,valid', [
    ('AACAA', 102, 'C', 'T', 'AATAA', True),
    ('AACAA', 102, 'G', 'T', 'AATAA', False),  # unexpected reference nucleotide
    ('AACAA', 1000, 'C', 'T', 'AATAA', False)  # variant out of range
])
def test_pam_protected_sequence_from_reference_sequence(seq, pos, ref, alt, ppseq, valid):
    chromosome = 'X'
    start = 100
    end = start + len(seq) - 1

    variant = PamVariant(GenomicPosition(chromosome, pos), ref, alt, 'sgrna-1')

    gr = GenomicRange(chromosome, start, end, '+')
    ref_seq = ReferenceSequence(seq, gr)
    with pytest.raises(Exception) if not valid else nullcontext():
        pam_ref_seq = PamBgAltSeqBuilder.from_ref_seq(
            ref_seq, [], [variant])
        alt_seq = pam_ref_seq.get_pam_seq(ref_check=True)

    if valid:
        assert pam_ref_seq.ref_seq == ref_seq.ref_seq
        assert alt_seq == ppseq


def test_pam_bg_alt_seq_builder_pam_seq():
    b = PamBgAltSeqBuilder.from_ref_seq(
        get_ref_seq('AAA', pos=100), [], [
            PamVariant(GenomicPosition('X', 101), 'A', 'C', 'sgrna-id-a')
        ])

    exp = 'ACA'
    assert b.get_pam_seq(extend=False) == exp
    assert b.pam_seq == exp
