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
from valiant.models.targeton import CDSTargeton
from .utils import get_aux_tables, seq2triplets

CONST_CODON_METHODS = {
    'A': 'get_ala_mutations',
    'STOP': 'get_stop_mutations'
}

aux = get_aux_tables()
sub_var_type = np.int8(VariantType.SUBSTITUTION.value)


@pytest.mark.parametrize('aa,strand,seq,exp_mseq', [
    ('A', '+', 'GCAAAATTT', ['GCCAAATTT', 'GCAGCCTTT', 'GCAAAAGCC']),
    ('A', '+', 'GCCAAATTT', ['GCCGCCTTT', 'GCCAAAGCC']),
    ('A', '-', 'TGCAAATTT', ['GGCAAATTT', 'TGCGGCTTT', 'TGCAAAGGC']),
    ('A', '+', 'GATTTT', ['GCCTTT', 'GATGCC']),
    ('STOP', '+', 'TAACCCGGG', ['TGACCCGGG', 'TAATGAGGG', 'TAACCCTGA']),
    ('STOP', '-', 'TTACCCGGG', ['TCACCCGGG', 'TTATCAGGG', 'TTACCCTCA'])
])
def test_get_snvre_aa_mutations(aa, strand, seq, exp_mseq):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, strand)
    t = CDSTargeton.from_pam_seq(PamProtectedReferenceSequence(seq, gr, seq), '', '')

    # Generate codon substitution mutations
    mc = getattr(t, CONST_CODON_METHODS[aa])(aux_tables=aux)

    # Check results
    assert mc.df.alt_aa.unique()[0] == aa
    assert np.array_equal(mc.df.mseq.to_numpy(), np.array(exp_mseq))
    assert mc.df.var_type.unique()[0] == sub_var_type
    assert np.array_equal(
        mc.df.ref.astype('string').to_numpy(),
        np.array([
            triplet
            for i, triplet in enumerate(seq2triplets(seq))
            if 3 * i in mc.df.mut_position.unique()
        ]))
