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

from hashlib import md5
import os
import pathlib
import pandas as pd
from valiant.models.base import GenomicRange
from valiant.models.codon_table import CodonTable
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.sequences import ReferenceSequence
from valiant.models.snv_table import AuxiliaryTables
from valiant.models.targeton import Targeton
from .constants import DUMMY_PAM_PROTECTION_NT, CODON_TABLE_FP, FRAMES, STRANDS

cds_extension_slices = {
    (a, b): slice(a, -b if b > 0 else None)
    for a in FRAMES
    for b in FRAMES
}


def get_data_file_path(fp):
    return os.path.join(pathlib.Path(__file__).parent.absolute(), 'data', fp)


def get_dummy_pam_protected(seq):
    return DUMMY_PAM_PROTECTION_NT * len(seq)


def get_pam_protected_sequence(seq, pam_protection, chromosome='X', strand='+', pos=1):
    gr = GenomicRange(chromosome, pos, len(seq), strand)
    ref_seq = ReferenceSequence(seq, gr)
    return PamProtectedReferenceSequence.from_reference_sequence(
        ref_seq, get_dummy_pam_protected(seq) if pam_protection else seq)


def get_targeton(seq, pam_protection, chromosome='X', strand='+', pos=1):
    pam_ref_seq = get_pam_protected_sequence(
        seq, pam_protection, chromosome=chromosome, strand=strand, pos=pos)
    return Targeton.from_pam_seq(pam_ref_seq)


def load_codon_table():
    return CodonTable.load(get_data_file_path(CODON_TABLE_FP))


def get_aux_tables(strands=None):
    strands_ = strands if strands is not None else STRANDS
    return AuxiliaryTables(load_codon_table(), strands_, True, True, True)


def load_csv(fp):
    return pd.read_csv(get_data_file_path(fp))


def md5sum(s):
    return md5(s.encode('utf-8')).hexdigest()


def trim_cds_extensions(cds_seq, cds_prefix_length, cds_suffix_length):
    assert len(cds_seq) % 3 == 0
    return cds_seq[cds_extension_slices[cds_prefix_length, cds_suffix_length]]


def seq2triplets(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]
