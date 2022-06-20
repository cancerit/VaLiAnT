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
from valiant.models.base import GenomicRange
from valiant.models.refseq_repository import ReferenceSequenceRepository, _normalise_ref_seq


@pytest.mark.parametrize('seq,norm_seq', [
    ('AACGT', 'AACGT'),
    ('aaCgT', 'AACGT'),  # soft-masked sequence
    ('XXACGT', None)     # invalid symbols
])
def test_normalise_ref_seq(seq, norm_seq):
    if norm_seq is not None:
        assert _normalise_ref_seq(seq) == norm_seq
    else:
        with pytest.raises(ValueError):
            _normalise_ref_seq(seq)


def test_reference_sequence_repository_register_genomic_ranges():
    rsr = ReferenceSequenceRepository()
    grs = [
        GenomicRange('X', 1, 100, '+'),
        GenomicRange('X', 500, 800, '+')
    ]
    rsr.register_genomic_ranges(grs)
    assert rsr.region_count == len(grs)


def test_reference_sequence_repository_register_get_sequence():
    chromosome = 'X'
    start = 1
    end = 100
    seq = 'AACCGGTT'
    gr = GenomicRange(chromosome, start, end, '+')

    # Initialise repository
    rsr = ReferenceSequenceRepository()
    rsr.register_genomic_range(gr)

    # Register sequence
    rsr.register_sequence(chromosome, start, end, seq)

    # Retrieve sequence
    assert rsr.get_sequence(chromosome, start, end) == seq
    assert rsr.get_genomic_range_sequence(gr) == seq


def test_reference_sequence_repository_get_genomic_range_subsequence():
    genomic_range = GenomicRange('X', 1023, 1032, '+')
    a = 'AAAAA'
    b = 'GGGGG'
    sequence = a + b
    rsr = ReferenceSequenceRepository()
    rsr.register_genomic_range(genomic_range)
    rsr._sequences[genomic_range.chromosome][(genomic_range.start, genomic_range.end)] = sequence

    assert rsr.get_genomic_range_subsequence(genomic_range, 1023, 1027) == a
    assert rsr.get_genomic_range_subsequence(genomic_range, 1028, 1032) == b
    assert a + b == sequence
