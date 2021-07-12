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

import logging
import sys
from valiant.models.base import GenomicRange
from valiant.models.refseq_repository import ReferenceSequenceRepository, fetch_reference_sequences, get_fasta_file


def test_reference_sequence_retrieval(ref_fp):
    ref = ReferenceSequenceRepository()
    fasta_file = get_fasta_file(ref_fp)

    chromosome = 'X'
    start = 41341615
    end = 41341625

    r = (chromosome, start, end)

    # Register region
    ref.register_region(*r)

    # Check no sequence has been retrieved yet
    seq = ref.get_sequence(*r)
    assert seq is None

    # Retrieve sequences for all regions
    ref.fetch_sequences(fasta_file)
    fasta_file.close()

    # Check the sequence has been retrieved
    seq = ref.get_sequence(*r)
    print(seq)
    print(len(seq))
    assert len(seq) == end - start + 1


def test_reference_subsequence(ref_fp):
    ref = ReferenceSequenceRepository()
    fasta_file = get_fasta_file(ref_fp)

    chromosome = 'X'
    strand = '+'
    start = 41341615
    end = 41341635

    ref_range = GenomicRange(chromosome, start, end, strand)

    ref.register_genomic_range(ref_range)
    ref.fetch_sequences(fasta_file)
    seq = ref.get_genomic_range_sequence(ref_range)
    offset = 10
    pre = ref.get_genomic_range_subsequence(ref_range, start, start + offset)
    print(seq)
    print(pre)
    print(len(pre))
    assert pre == seq[:offset]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the path to the reference genome FASTA file!")
        sys.exit(1)

    logging.basicConfig(level=logging.DEBUG)

    ref_fp = sys.argv[1]
    test_reference_sequence_retrieval(ref_fp)
    test_reference_subsequence(ref_fp)
