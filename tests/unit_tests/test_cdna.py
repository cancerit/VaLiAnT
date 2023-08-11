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
from functools import lru_cache
from valiant.enums import TargetonMutator
from valiant.models.base import PositionRange, StrandedPositionRange, GenomicRange
from valiant.models.cdna import CDNA, AnnotatedCDNA
from valiant.models.cdna_targeton_configs import CDNATargetonConfig, CDNATargetonConfigCollection
from valiant.models.sequences import Sequence
from valiant.models.sequence_info import SequenceInfo
from .utils import get_data_file_path


cdna_targeton_fp = get_data_file_path("cdna_targeton.tsv")

SEQ = "ATGCCCTGA"


@lru_cache(maxsize=4)
def get_cdna_annot(seq, cds_range=None):
    return AnnotatedCDNA(
        Sequence(seq),
        SequenceInfo.empty(),
        cds_range or StrandedPositionRange(1, len(seq), '+'))


def test_cdna_targeton_collection_load():
    cdna_tcfgs = CDNATargetonConfigCollection.load(cdna_targeton_fp)
    assert len(cdna_tcfgs) == 1

    tcfg = cdna_tcfgs.cts[0]
    assert tcfg.seq_id == 'cdna1'
    assert tcfg.targeton_range.start == 10
    assert tcfg.targeton_range.end == 80
    assert tcfg.r2_range.start == 20
    assert tcfg.r2_range.end == 70
    assert tcfg.mutators == {TargetonMutator.SNV, TargetonMutator.DEL1}


@pytest.mark.parametrize('start,end,exp_subseq,exp_ext_5p,exp_ext_3p', [
    (3, 6, 'GCCC', 'AT', ''),
    (4, 6, 'CCC', '', ''),
    (6, 7, 'CT', 'CC', 'GA')
])
def test_annotated_cdna_get_extended_subsequence(start, end, exp_subseq, exp_ext_5p, exp_ext_3p):
    cdna = get_cdna_annot(SEQ)
    subseq, ext_5p, ext_3p = cdna.get_extended_subsequence(PositionRange(start, end))
    assert subseq.sequence == exp_subseq
    assert ext_5p == exp_ext_5p
    assert ext_3p == exp_ext_3p
