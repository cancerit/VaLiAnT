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

from dataclasses import dataclass
import pytest
from valiant.enums import TargetonMutator
from valiant.models.base import GenomicRange, TranscriptInfo
from valiant.models.codon_table import CodonTable
from valiant.models.oligo_template import OligoTemplate, TargetonOligoSegment
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.sequences import Sequence, ReferenceSequence
from valiant.models.targeton import Targeton
from .constants import CODON_TABLE_FP, DUMMY_PAM_PROTECTION_NT
from .utils import get_data_file_path, get_targeton

TRANSCRIPT_INFO = TranscriptInfo('GENE_ID', 'TRANSCRIPT_ID')


@pytest.mark.parametrize('targetons', [
    ['GGGACGTACGTTTTTACGTGGG', 'ATATATATATATATATATATATA', 'GGGGTATGCGGATTGACTGCTAGT']
])
@pytest.mark.parametrize('pam_protection', [True, False])
@pytest.mark.parametrize('mutator', [TargetonMutator.DEL1, TargetonMutator.SNV])
def test_oligo_compute_mutations(targetons, mutator, pam_protection):

    def get_segment(seq):
        return TargetonOligoSegment(
            get_targeton(seq, pam_protection), {mutator})

    ct = CodonTable.load(get_data_file_path(CODON_TABLE_FP))
    ref_seq = ''.join(targetons)
    pam_ref_seq = PamProtectedReferenceSequence(
        ref_seq,
        GenomicRange('X', 1, sum(len(seq) for seq in targetons), '+'),
        ref_seq)

    adaptor_5 = 'AAAAAA'
    adaptor_3 = 'AAAAAA'
    segments = list(map(get_segment, targetons))

    ot = OligoTemplate(TRANSCRIPT_INFO, pam_ref_seq, set(), set(), adaptor_5, adaptor_3, segments)
    for _, target_segment in ot.target_segments:
        mutation_collections = target_segment.compute_mutations(ct)
        mutation_collection = mutation_collections[mutator]
        if pam_protection:
            if mutator == TargetonMutator.DEL1:
                assert all(
                    set(m.sequence) == {DUMMY_PAM_PROTECTION_NT}
                    for m in mutation_collection.mutations
                )
            elif mutator == TargetonMutator.SNV:
                assert all(
                    set(m.sequence) - {m.new} == {DUMMY_PAM_PROTECTION_NT}
                    for m in mutation_collection.mutations
                )
