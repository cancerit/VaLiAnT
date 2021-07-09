########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
#
# The Software remains the property of Genome Research Ltd ("GRL").
#
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
#
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
#
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
#
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
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
