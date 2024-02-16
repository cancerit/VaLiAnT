########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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
from valiant.cds_seq import CdsSeq
from valiant.exon import Exon
from valiant.seq import Seq
from valiant.strings.dna_str import DnaStr
from valiant.strings.strand import Strand
from valiant.transcript import Transcript
from valiant.transcript_info import TranscriptInfo
from valiant.uint_range import UIntRangeSortedList


seq = Seq(100, DnaStr('AAACCTTCGGGAATACCC'))
exons = [
    Exon(100, 104, 0, 0),
    Exon(107, 112, 1, 1),
    Exon(114, 117, 2, 1)
]


def get_transcript(strand, exons):
    return Transcript(
        TranscriptInfo('X', Strand(strand), None, None),
        UIntRangeSortedList(exons))


transcript_plus = get_transcript('+', exons)

codons = [
    CdsSeq(100, DnaStr('AAA')),
    CdsSeq(103, DnaStr('CC'), cds_suffix=DnaStr('C')),
    CdsSeq(108, DnaStr('GGG')),
    CdsSeq(111, DnaStr('AA'), cds_suffix=DnaStr('A')),
    CdsSeq(115, DnaStr('CCC'))
]


@pytest.mark.parametrize('pos,exp', [
    (50, None),   # out of bounds
    (105, None),  # intronic
    (200, None),  # out of bounds
    (100, codons[0]),
    (102, codons[0]),
    (103, codons[1]),
    (107, codons[1])
])
def test_transcript_get_codon_at(pos, exp):
    obs = transcript_plus.get_codon_at(seq, pos)
    if exp is None:
        assert obs is None
    else:
        assert obs is not None
        assert obs == exp or obs.ext == exp.ext
