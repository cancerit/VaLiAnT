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

from contextlib import nullcontext
import pytest
from valiant.errors import SequenceNotFound
from valiant.loaders.fasta import load_from_multi_fasta
from .constants import MULTI_FASTA_FP
from .utils import get_data_file_path


SEQ_ID_1 = 'seq001'
SEQ_ID_2 = 'seq002'

multi_fasta_fp = get_data_file_path(MULTI_FASTA_FP)


@pytest.mark.parametrize('seq_ids,all_exist', [
    ([SEQ_ID_1, SEQ_ID_2], True),
    ([SEQ_ID_1, 'missing-id'], False)
])
def test_load_from_multi_fasta(seq_ids, all_exist):
    with pytest.raises(SequenceNotFound) if not all_exist else nullcontext():
        id_seqs = load_from_multi_fasta(multi_fasta_fp, seq_ids)

    if all_exist:
        assert set(seq_ids) == set(id_seqs.keys())
