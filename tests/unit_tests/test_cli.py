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
from valiant.cli import _load_codon_table, _load_gff_file, _load_pam_protection_vcf
from valiant.models.codon_table import CodonTable
from valiant.models.exon import AnnotationRepository
from .constants import CODON_TABLE_FP, GTF_SINGLE, GTF_MULTI, PAM_VCF_FP
from .utils import get_data_file_path


def test_load_codon_table_default():

    # Load default codon table
    default_codon_table = _load_codon_table(None)
    assert isinstance(default_codon_table, CodonTable)
    aas = set(default_codon_table.amino_acid_symbols)
    assert 'STOP' in aas
    assert len(aas) == 21


def test_load_codon_table_invalid():

    # Load invalid codon table
    with pytest.raises(SystemExit):
        _load_codon_table(get_data_file_path(GTF_SINGLE))


def test_load_codon_table_valid():

    # Load valid codon table
    codon_table = _load_codon_table(get_data_file_path(CODON_TABLE_FP))

    # Check codon table
    assert isinstance(codon_table, CodonTable)


def test_load_gff_file_none():

    # Load no GTF
    assert _load_gff_file(None) is None


def test_load_gff_file_invalid():
    with pytest.raises(SystemExit):

        # Load invalid GTF
        _load_gff_file(get_data_file_path(GTF_MULTI))


def test_load_gff_file_valid():

    # Load valid GTF
    annotation = _load_gff_file(get_data_file_path(GTF_SINGLE))

    # Check annotation repository
    assert isinstance(annotation, AnnotationRepository)
    assert annotation.cds is not None


def test_load_pam_protection_vcf_none():

    # Load no PAM VCF file
    vr = _load_pam_protection_vcf(set(), None)
    assert vr.sgrna_ids == set()


def test_load_pam_protection_vcf_invalid():

    # Load invalid PAM VCF file
    with pytest.raises(SystemExit):
        _load_codon_table(get_data_file_path(GTF_SINGLE))


def test_load_pam_protection_vcf_valid():
    sgrna_ids = {'sgRNA_1', 'sgRNA_2', 'sgRNA_3', 'sgRNA_4'}
    pam_fp = get_data_file_path(PAM_VCF_FP)

    # Load valid PAM VCF file
    vr = _load_pam_protection_vcf(sgrna_ids, pam_fp)

    # Check PAM variant repository
    assert vr.sgrna_ids == sgrna_ids
    assert len(vr._variants) == 4
    assert len(vr._ranges) == 5
