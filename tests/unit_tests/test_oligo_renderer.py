########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
