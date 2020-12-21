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

import logging
import sys
from valiant.models.base import GenomicRange
from valiant.models.refseq_repository import ReferenceSequenceRepository, fetch_reference_sequences, get_fasta_file


def test_reference_sequence_retrieval(ref_fp):
    ref = ReferenceSequenceRepository()
    fasta_file = get_fasta_file(ref_fp)

    chromosome = 'X'
    start = 41341615
    end = 41341625

    r = (chromosome, start, end)

    # Register region
    ref.register_region(*r)

    # Check no sequence has been retrieved yet
    seq = ref.get_sequence(*r)
    assert seq is None

    # Retrieve sequences for all regions
    ref.fetch_sequences(fasta_file)
    fasta_file.close()

    # Check the sequence has been retrieved
    seq = ref.get_sequence(*r)
    print(seq)
    print(len(seq))
    assert len(seq) == end - start + 1


def test_reference_subsequence(ref_fp):
    ref = ReferenceSequenceRepository()
    fasta_file = get_fasta_file(ref_fp)

    chromosome = 'X'
    strand = '+'
    start = 41341615
    end = 41341635

    ref_range = GenomicRange(chromosome, start, end, strand)

    ref.register_genomic_range(ref_range)
    ref.fetch_sequences(fasta_file)
    seq = ref.get_genomic_range_sequence(ref_range)
    offset = 10
    pre = ref.get_genomic_range_subsequence(ref_range, start, start + offset)
    print(seq)
    print(pre)
    print(len(pre))
    assert pre == seq[:offset]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the path to the reference genome FASTA file!")
        sys.exit(1)

    logging.basicConfig(level=logging.DEBUG)

    ref_fp = sys.argv[1]
    test_reference_sequence_retrieval(ref_fp)
    test_reference_subsequence(ref_fp)
