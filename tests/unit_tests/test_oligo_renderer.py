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
import pandas as pd
import pytest
from valiant.enums import VariantType, TargetonMutator
from valiant.models.base import GenomicRange
from valiant.models.options import Options
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.oligo_renderer import BaseOligoRenderer
from .utils import get_pam_protected_sequence


SEQ = 'AAAACGTACGTACGT'
GENE_ID = 'G1'
TRANSCRIPT_ID = 'T1'
META_COLUMNS = {
    'gene_id',
    'transcript_id',
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'ref_seq',
    'pam_seq',
    'revc'
}


def test_base_oligo_renderer_init():
    pam_seq = get_pam_protected_sequence(SEQ, None)

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, 'G1', 'T1', '', '')

    # Check renderer properties
    assert renderer.chromosome == pam_seq.genomic_range.chromosome
    assert renderer.start == pam_seq.genomic_range.start
    assert renderer.strand == pam_seq.genomic_range.strand


def test_base_oligo_renderer_constant_fields():
    gene_id = GENE_ID
    transcript_id = TRANSCRIPT_ID
    pam_seq = get_pam_protected_sequence(SEQ, None)

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, gene_id, transcript_id, '', '')

    # Check constant fields
    const_fields = renderer._constant_fields
    assert const_fields[0] == ('gene_id', gene_id)
    assert const_fields[1] == ('transcript_id', transcript_id)
    assert const_fields[2] == ('ref_chr', renderer.chromosome)
    assert const_fields[3] == ('ref_strand', renderer.ref_seq.genomic_range.strand)
    assert const_fields[4] == ('ref_seq', renderer.ref_seq.sequence)
    assert const_fields[5] == ('pam_seq', renderer.ref_seq.pam_protected_sequence)


@pytest.mark.parametrize('seq,a5,a3,exp,exp_rc', [
    ('AACCGGCC', 'TTT', 'TTT', 'TTTAACCGGCCTTT', 'TTTGGCCGGTTTTT')
])
def test_base_oligo_renderer_render_mutated_sequence(seq, a5, a3, exp, exp_rc):
    pam_seq = get_pam_protected_sequence(SEQ, None)

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, GENE_ID, TRANSCRIPT_ID, a5, a3)

    # Check rendered oligonucleotide sequence
    assert renderer._render_mutated_sequence(seq) == exp

    # Check rendered oligonucleotide sequence (reverse complement)
    assert renderer._render_mutated_sequence_rc(seq) == exp_rc


@pytest.mark.parametrize('strand,rc,exp,valid', [
    ('+', False, '_render_mutated_sequence', True),
    ('+', True, None, False),
    ('-', False, '_render_mutated_sequence', True),
    ('-', True, '_render_mutated_sequence_rc', True)
])
def test_base_oligo_renderer_get_renderer(strand, rc, exp, valid):
    pam_seq = get_pam_protected_sequence(SEQ, None, strand=strand)

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, GENE_ID, TRANSCRIPT_ID, '', '')

    with pytest.raises(ValueError) if not valid else nullcontext():
        assert renderer._get_renderer(rc).__name__ == exp


@pytest.mark.parametrize('var_type,source,start,ref,alt,exp', [
    (VariantType.DELETION, TargetonMutator.DEL1.value, 100, 'A', None, 'T1.G1_X:100_1del'),
    (VariantType.DELETION, TargetonMutator.DEL2_0.value, 100, 'AA', None, 'T1.G1_X:100_101_2del0'),
    (VariantType.DELETION, TargetonMutator.DEL2_1.value, 100, 'AA', None, 'T1.G1_X:100_101_2del1'),
    (VariantType.INSERTION, 'vcf', 100, None, 'ACGT', 'T1.G1_X:100_ACGT_vcf'),
    (VariantType.SUBSTITUTION, TargetonMutator.SNV.value, 100, 'AAA', 'CGT', 'T1.G1_X:100_102_AAA>CGT_snv'),
    (VariantType.SUBSTITUTION, TargetonMutator.SNV.value, 100, 'A', 'C', 'T1.G1_X:100_A>C_snv'),
    (VariantType.SUBSTITUTION, 'vcf', 100, 'A', 'C', 'T1.G1_X:100_A>C_vcf'),
])
def test_base_oligo_renderer_get_oligo_name(var_type, source, start, ref, alt, exp):
    pam_seq = get_pam_protected_sequence(SEQ, None)

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, GENE_ID, TRANSCRIPT_ID, '', '')

    assert renderer._get_oligo_name(var_type.value, source, start, ref, alt) == exp


@pytest.mark.parametrize('rc', [True, False])
def test_base_oligo_renderer_get_metadata_table(rc):
    pam_seq = get_pam_protected_sequence(SEQ, None, strand='-')

    df = pd.DataFrame.from_records([
        ('OLIGO_1', 100, 'A', 'C', 'snv', 'AAACGGG')
    ], columns=['oligo_name', 'mut_position', 'ref', 'new', 'mutator', 'mseq'])

    # Initialise renderer
    renderer = BaseOligoRenderer(pam_seq, GENE_ID, TRANSCRIPT_ID, '', '')

    meta = renderer.get_metadata_table(df, Options(rc, 300))
    meta_cols = set(meta.columns)
    assert META_COLUMNS <= meta_cols

    name_rc_mask = meta.oligo_name.str.endswith('_rc')
    revc_set = set(meta.revc.unique())

    if rc:
        assert name_rc_mask.all()
        assert revc_set == {1}
    else:
        assert not name_rc_mask.any()
        assert revc_set == {0}
