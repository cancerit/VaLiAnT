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
from valiant.string_mutators import (
    replace_single_nucleotides,
    delete_non_overlapping_2,
    delete_non_overlapping_3_offset
)
from .constants import DNA_ALPHABET, DNA_ALPHABET_SET


@pytest.mark.parametrize('seq,exp_res', [
    ('CC', {
        (0, 'C', 'A', 'AC'),
        (0, 'C', 'G', 'GC'),
        (0, 'C', 'T', 'TC'),
        (1, 'C', 'A', 'CA'),
        (1, 'C', 'G', 'CG'),
        (1, 'C', 'T', 'CT')
    }),
    ('AAA', {
        (0, 'A', 'C', 'CAA'),
        (0, 'A', 'G', 'GAA'),
        (0, 'A', 'T', 'TAA'),
        (1, 'A', 'C', 'ACA'),
        (1, 'A', 'G', 'AGA'),
        (1, 'A', 'T', 'ATA'),
        (2, 'A', 'C', 'AAC'),
        (2, 'A', 'G', 'AAG'),
        (2, 'A', 'T', 'AAT')
    })
])
def test_replace_single_nucleotides(seq, exp_res):
    seq_n = len(seq)
    exp_mutation_n = 3 * seq_n

    # Generate mutations
    res = replace_single_nucleotides(seq, list(DNA_ALPHABET_SET))

    # Validate resulting number of mutations
    assert len(res) == exp_mutation_n

    # Validate mutation items format
    for i, ref_nt, new_nt, mseq in res:
        assert isinstance(i, int)
        assert ref_nt in DNA_ALPHABET_SET
        assert new_nt in DNA_ALPHABET_SET
        assert len(mseq) == seq_n
        assert set(mseq) <= DNA_ALPHABET_SET
        assert mseq != seq

    # Validate full set of mutations
    assert set(res) == exp_res


@pytest.mark.parametrize('seq,offset,exp_res', [
    ('AAACCC', 0, {
        (0, 'AA', 'ACCC'),
        (2, 'AC', 'AACC'),
        (4, 'CC', 'AAAC')
    }),
    ('AAACCCT', 0, {
        (0, 'AA', 'ACCCT'),
        (2, 'AC', 'AACCT'),
        (4, 'CC', 'AAACT')
    }),
    ('TAAACCC', 1, {
        (1, 'AA', 'TACCC'),
        (3, 'AC', 'TAACC'),
        (5, 'CC', 'TAAAC')
    }),
    ('AAACCC', 1, {
        (1, 'AA', 'ACCC'),
        (3, 'CC', 'AAAC')
    }),
    ('A', 0, set()),
    ('A', 1, set()),
    ('AA', 0, {
        (0, 'AA', '')
    })
])
def test_delete_non_overlapping_2(seq, offset, exp_res):
    res = delete_non_overlapping_2(seq, offset)
    assert set(res) == exp_res


@pytest.mark.parametrize('seq, const_prefix_len, const_suffix_len, exp_res', [
    ('AAACCC', 1, 2, {
        (1, 'AAC', 'ACC'),
    }),
    ('AAACCCGGGTTT', 1, 2, {
        (1, 'AAC', 'ACCGGGTTT'),
        (4, 'CCG', 'AAACGGTTT'),
        (7, 'GGT', 'AAACCCGTT')
    }),
    ('AAACCCGGGT', 1, 0, {
        (1, 'AAC', 'ACCGGGT'),
        (4, 'CCG', 'AAACGGT'),
        (7, 'GGT', 'AAACCCG')
    })
])
def test_delete_non_overlapping_3_offset(seq, const_prefix_len, const_suffix_len, exp_res):

    seq_n = len(seq)
    assert (seq_n - const_prefix_len - const_suffix_len) % 3 == 0

    res = delete_non_overlapping_3_offset(seq, const_prefix_len, const_suffix_len)
    for i, ref, mseq in res:
        print(i, ref, mseq)
        assert isinstance(i, int)
        assert len(ref) == 3
        assert set(ref) <= DNA_ALPHABET_SET
        assert set(mseq) <= DNA_ALPHABET_SET
        assert mseq != seq

    # Validate full set of mutations
    assert set(res) == exp_res
