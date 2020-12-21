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
from valiant.string_mutators import (
    replace_single_nucleotides,
    delete_non_overlapping_2,
    delete_non_overlapping_3_offset
)
from .constants import DNA_ALPHABET, DNA_ALPHABET_SET


@pytest.mark.parametrize('seq,exp_res', [
    ('CC', {
        (0, 'C', 'A', 'AC'),
        (0, 'C', 'G', 'GC'),
        (0, 'C', 'T', 'TC'),
        (1, 'C', 'A', 'CA'),
        (1, 'C', 'G', 'CG'),
        (1, 'C', 'T', 'CT')
    }),
    ('AAA', {
        (0, 'A', 'C', 'CAA'),
        (0, 'A', 'G', 'GAA'),
        (0, 'A', 'T', 'TAA'),
        (1, 'A', 'C', 'ACA'),
        (1, 'A', 'G', 'AGA'),
        (1, 'A', 'T', 'ATA'),
        (2, 'A', 'C', 'AAC'),
        (2, 'A', 'G', 'AAG'),
        (2, 'A', 'T', 'AAT')
    })
])
def test_replace_single_nucleotides(seq, exp_res):
    seq_n = len(seq)
    exp_mutation_n = 3 * seq_n

    # Generate mutations
    res = replace_single_nucleotides(seq, list(DNA_ALPHABET_SET))

    # Validate resulting number of mutations
    assert len(res) == exp_mutation_n

    # Validate mutation items format
    for i, ref_nt, new_nt, mseq in res:
        assert isinstance(i, int)
        assert ref_nt in DNA_ALPHABET_SET
        assert new_nt in DNA_ALPHABET_SET
        assert len(mseq) == seq_n
        assert set(mseq) <= DNA_ALPHABET_SET
        assert mseq != seq

    # Validate full set of mutations
    assert set(res) == exp_res


@pytest.mark.parametrize('seq,offset,exp_res', [
    ('AAACCC', 0, {
        (0, 'AA', 'ACCC'),
        (2, 'AC', 'AACC'),
        (4, 'CC', 'AAAC')
    }),
    ('AAACCCT', 0, {
        (0, 'AA', 'ACCCT'),
        (2, 'AC', 'AACCT'),
        (4, 'CC', 'AAACT')
    }),
    ('TAAACCC', 1, {
        (1, 'AA', 'TACCC'),
        (3, 'AC', 'TAACC'),
        (5, 'CC', 'TAAAC')
    }),
    ('AAACCC', 1, {
        (1, 'AA', 'ACCC'),
        (3, 'CC', 'AAAC')
    }),
    ('A', 0, set()),
    ('A', 1, set()),
    ('AA', 0, {
        (0, 'AA', '')
    })
])
def test_delete_non_overlapping_2(seq, offset, exp_res):
    res = delete_non_overlapping_2(seq, offset)
    assert set(res) == exp_res


@pytest.mark.parametrize('seq, const_prefix_len, const_suffix_len, exp_res', [
    ('AAACCC', 1, 2, {
        (1, 'AAC', 'ACC'),
    }),
    ('AAACCCGGGTTT', 1, 2, {
        (1, 'AAC', 'ACCGGGTTT'),
        (4, 'CCG', 'AAACGGTTT'),
        (7, 'GGT', 'AAACCCGTT')
    }),
    ('AAACCCGGGT', 1, 0, {
        (1, 'AAC', 'ACCGGGT'),
        (4, 'CCG', 'AAACGGT'),
        (7, 'GGT', 'AAACCCG')
    })
])
def test_delete_non_overlapping_3_offset(seq, const_prefix_len, const_suffix_len, exp_res):

    seq_n = len(seq)
    assert (seq_n - const_prefix_len - const_suffix_len) % 3 == 0

    res = delete_non_overlapping_3_offset(seq, const_prefix_len, const_suffix_len)
    for i, ref, mseq in res:
        print(i, ref, mseq)
        assert isinstance(i, int)
        assert len(ref) == 3
        assert set(ref) <= DNA_ALPHABET_SET
        assert set(mseq) <= DNA_ALPHABET_SET
        assert mseq != seq

    # Validate full set of mutations
    assert set(res) == exp_res
