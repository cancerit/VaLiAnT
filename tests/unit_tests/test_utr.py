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

import pandas as pd
from pyranges import from_string
import pytest
from valiant.models.base import GenomicRange, TranscriptInfo
from valiant.models.exon import UTRRepository

UTR_RANGES = from_string("""
Chromosome	Start	End	Strand	gene_id	transcript_id
chrX	1000	1300	+	G1	T1
chrX	5000	5300	-	G2	T2
""")

REF_RANGES = from_string("""
Chromosome	Start	End	Strand
chrX	799	1200	+
chrX	4899	5400	+
""")


def test_utr_repository_init():
    utr = UTRRepository(UTR_RANGES)
    assert len(utr.utr_ranges) == 2
    assert utr.utr_ranges.chromosomes == ['chrX']


def test_utr_repository_get_transcript_infos():
    utr = UTRRepository(UTR_RANGES)
    d = utr.get_transcript_infos(REF_RANGES)
    assert len(d) == 1
    gr = GenomicRange('chrX', 800, 1200, '+')
    ti = TranscriptInfo('G1', 'T1')
    assert d[gr] == ti
