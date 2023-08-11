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
from valiant.models.base import GenomicRange
from valiant.models.sequences import Sequence, ReferenceSequence


@pytest.mark.parametrize('seq,valid', [
    ('AAA', True),
    ('Aaa', False),
    ('', False)
])
def test_sequence_init(seq, valid):
    with pytest.raises(ValueError) if not valid else nullcontext():

        # Initialise sequence
        sequence = Sequence(seq)

    if valid:

        # Check sequence length
        assert len(seq) == len(sequence)


@pytest.mark.parametrize('seq,start,end,valid', [
    ('AAA', 10, 12, True),
    ('AAA', 10, 13, False)
])
def test_reference_sequence_init(seq, start, end, valid):
    gr = GenomicRange('X', start, end, '+')

    with pytest.raises(ValueError) if not valid else nullcontext():

        # Initialise reference sequence
        ReferenceSequence(seq, gr)


def test_reference_sequence_get_subsequence():
    seq = 'AAACCC'
    gr = GenomicRange('X', 100, 105, '+')
    gr_sub = GenomicRange('X', 100, 102, '+')

    # Initialise reference sequence
    ref_seq = ReferenceSequence(seq, gr)

    # Extract reference subsequence
    sub_ref_seq = ref_seq.get_subsequence(gr_sub)

    # Check reference subsequence
    assert sub_ref_seq.genomic_range == gr_sub
    assert sub_ref_seq.ref_seq == 'AAA'
