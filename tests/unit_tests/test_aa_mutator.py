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
from valiant.models.codon_table import STOP_CODE
from valiant.snv import build_aa_sub_table, get_fast_aa_subs
from .utils import load_codon_table, get_aux_tables

codon_table = load_codon_table()


@pytest.mark.parametrize('strand', ['+', '-'])
def test_build_aa_sub_table(strand):
    t = build_aa_sub_table(codon_table, strand)
    assert t.shape[0] == 1219
    t.ref_aa.cat.remove_unused_categories(inplace=True)
    t.alt_aa.cat.remove_unused_categories(inplace=True)
    assert STOP_CODE in t.ref_aa.cat.categories
    assert STOP_CODE not in t.alt_aa.cat.categories


aux = get_aux_tables()
all_aa_table = aux.all_aa_table

@pytest.mark.parametrize('seq', ['AAACCCTGA'])
@pytest.mark.parametrize('frame', [0, 1, 2])
@pytest.mark.parametrize('prefix', [0, 1, 2])
@pytest.mark.parametrize('suffix', [0, 1, 2])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_get_fast_aa_subs(seq, frame, prefix, suffix, strand):
    seq_ = ('T' * prefix) + seq + ('T' * suffix)
    df = get_fast_aa_subs(
        all_aa_table._strand_table[strand], frame, seq_)
    assert not df.ref.isna().any()
