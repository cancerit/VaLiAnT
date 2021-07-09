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

import numpy as np
import pytest
from valiant.enums import VariantType
from valiant.models.base import GenomicRange
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.targeton import Targeton, CDSTargeton

del_var_type = np.int8(VariantType.DELETION.value)

del_offset_method = {
    0: 'get_2del0_mutations',
    1: 'get_2del1_mutations'
}


@pytest.mark.parametrize('offset,seq,exp_pos,exp_ref,exp_mseq', [
    (0, 'AAATTT', [0, 2, 4], ['AA', 'AT', 'TT'], ['ATTT', 'AATT', 'AAAT']),
    (1, 'AAATTT', [1, 3], ['AA', 'TT'], ['ATTT', 'AAAT']),
    (0, 'AATTT', [0, 2], ['AA', 'TT'], ['TTT', 'AAT']),
    (1, 'AATTT', [1, 3], ['AT', 'TT'], ['ATT', 'AAT'])
])
@pytest.mark.parametrize('cds', [True, False])
def test_get_2del_mutations(offset, seq, exp_pos, exp_ref, exp_mseq, cds):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, '+')
    t = (
        CDSTargeton.from_pam_seq(PamProtectedReferenceSequence(seq, gr, seq), 'AA', 'A') if cds else
        Targeton.from_pam_seq(PamProtectedReferenceSequence(seq, gr, seq))
    )

    # Generate in-frame deletions
    mc = getattr(t, del_offset_method[offset])()

    # Check metadata table
    assert np.array_equal(mc.df.mut_position, np.array(exp_pos))
    assert np.array_equal(mc.df.ref.astype('string').to_numpy(), np.array(exp_ref))
    assert np.array_equal(mc.df.mseq.astype('string').to_numpy(), np.array(exp_mseq))
    assert mc.df.new.isna().all()
    assert mc.df.var_type.unique()[0] == del_var_type


@pytest.mark.parametrize('seq,pre,suf,strand,exp_pos,exp_ref,exp_mseq', [
    ('AAATTT', '', '', '+', [0, 3], ['AAA', 'TTT'], ['TTT', 'AAA']),
    ('GAAATTT', 'GG', '', '+', [1, 4], ['AAA', 'TTT'], ['GTTT', 'GAAA']),
    ('AAATTTGG', '', 'G', '+', [0, 3], ['AAA', 'TTT'], ['TTTGG', 'AAAGG']),
    ('GAAATTT', 'GG', '', '-', [1, 4], ['AAA', 'TTT'], ['GTTT', 'GAAA']),
    ('AAATTTGG', '', 'G', '-', [0, 3], ['AAA', 'TTT'], ['TTTGG', 'AAAGG'])
])
def test_get_inframe_mutations(seq, pre, suf, strand, exp_pos, exp_ref, exp_mseq):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, strand)
    t = CDSTargeton.from_pam_seq(PamProtectedReferenceSequence(seq, gr, seq), pre, suf)

    # Generate in-frame deletions
    mc = t.get_inframe_mutations()

    # Check metadata table
    assert np.array_equal(mc.df.mut_position, np.array(exp_pos))
    assert np.array_equal(mc.df.ref.astype('string').to_numpy(), np.array(exp_ref))
    assert np.array_equal(mc.df.mseq.astype('string').to_numpy(), np.array(exp_mseq))
    assert mc.df.new.isna().all()
    assert mc.df.var_type.unique()[0] == del_var_type
