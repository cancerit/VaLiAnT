########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from valiant.constants import META_PAM_MUT_SGRNA_ID, META_MUT_POSITION
from valiant.enums import MutationType, TargetonMutator
from valiant.models.annotated_sequence import CDSAnnotatedSequencePair
from valiant.models.base import GenomicPosition, GenomicRange
from valiant.models.pam_protection import PamVariant
from valiant.models.targeton import PamProtCDSTargeton
from .utils import get_aux_tables


aux = get_aux_tables()
codon_table = aux.codon_table


SGRNA_ID = 'some-id'
TEST_CHROMOSOME = 'X'


def test_pam_prot_targeton_from_reference_sequence():
    chromosome = TEST_CHROMOSOME
    PamProtCDSTargeton(CDSAnnotatedSequencePair(
        GenomicRange(chromosome, 500, 505, '+'),
        'ACGTCA',
        'AGGTCA',
        [PamVariant(GenomicPosition(chromosome, 501), 'C', 'G', SGRNA_ID)],
        'AT',
        'C'))


@pytest.mark.parametrize('alt,mut_type', [
    ('T', MutationType.SYNONYMOUS),
    ('G', MutationType.MISSENSE),
    ('A', MutationType.NONSENSE)
])
def test_pam_prot_targeton_get_pam_variant_annotations(alt, mut_type):
    # Codon table
    # TGA -> STOP
    # TGC -> C
    # TGT -> C
    # TGG -> W

    chromosome = TEST_CHROMOSOME
    PamProtCDSTargeton(CDSAnnotatedSequencePair(
        GenomicRange(chromosome, 500, 505, '+'),
        'GCGTGC',
        f"G{alt}GTGC",
        [PamVariant(GenomicPosition(chromosome, 501), 'C', alt, SGRNA_ID)],
        'T',
        'CA'))


@pytest.mark.parametrize('pam_variant_clash', [True, False])
def test_pam_prot_targeton_compute_mutations(pam_variant_clash):
    chromosome = TEST_CHROMOSOME
    ref_seq = 'ACGTCAGCG'
    genomic_range = GenomicRange(chromosome, 500, 500 + len(ref_seq) - 1, '+')
    mutator = TargetonMutator.SNV
    pam_variant_pos = 501 if pam_variant_clash else 507
    sgrna_id_1 = 'sgrna-a'
    sgrna_id_2 = 'sgrna-b'
    variants = [
        PamVariant(GenomicPosition(chromosome, 500), 'A', 'G', sgrna_id_1),
        PamVariant(GenomicPosition(chromosome, pam_variant_pos), 'C', 'G', sgrna_id_2)
    ]
    pam_seq = variants[0].mutate(ref_seq, genomic_range.start, ref_check=True)
    pam_seq = variants[1].mutate(pam_seq, genomic_range.start, ref_check=True)
    cds_prefix = 'A'
    frame = len(cds_prefix)

    # Generate targeton
    with pytest.raises(ValueError) if pam_variant_clash else nullcontext():
        targeton: PamProtCDSTargeton = PamProtCDSTargeton(
            CDSAnnotatedSequencePair(
                genomic_range, ref_seq, pam_seq, variants, cds_prefix, 'CG'))

    if not pam_variant_clash:

        # Generate mutations
        mutation_collections = targeton.compute_mutations(
            frozenset([mutator]), aux, sgrna_ids=frozenset([sgrna_id_1, sgrna_id_2]))

        assert len(mutation_collections) == 1

        # Verify the temporary field was dropped

        df = mutation_collections[mutator].df

        # Validate the sgRNA ID field
        df['codon_index'] = (df[META_MUT_POSITION] + frame) // 3
        mask = df['codon_index'].mod(2).ne(0)
        assert df.loc[mask, META_PAM_MUT_SGRNA_ID].unique() == ''
        assert df.loc[df.codon_index == 0, META_PAM_MUT_SGRNA_ID].unique() == sgrna_id_1
        assert df.loc[df.codon_index == 2, META_PAM_MUT_SGRNA_ID].unique() == sgrna_id_2


def test_pam_prot_cds_targeton_concat():
    def get_targeton(start, end, ref_seq, variants=None, alt_seq=None, cds_prefix='', cds_suffix=''):
        return PamProtCDSTargeton(
            CDSAnnotatedSequencePair(
                GenomicRange('X', start, end, '+'),
                ref_seq,
                alt_seq or ref_seq,
                variants or [],
                cds_prefix,
                cds_suffix))

    pam_variant = PamVariant(GenomicPosition('X', 47), 'A', 'T', 'sgrna-1')
    pam_variants = [pam_variant]

    t1 = get_targeton(45, 50, 'AATAAA', cds_prefix='T', cds_suffix='CC', variants=pam_variants)
    t2 = get_targeton(51, 56, 'CCCCCC', cds_prefix='G', cds_suffix='GG')
    t3 = get_targeton(57, 62, 'GGGGGG', cds_prefix='C', cds_suffix='AA')
    ct = PamProtCDSTargeton.concat([t1, t2, t3])

    assert ct.annotated_seq.cds_prefix == t1.annotated_seq.cds_prefix
    assert ct.annotated_seq.cds_suffix == t3.annotated_seq.cds_suffix
    assert ct.pos_range.start == t1.pos_range.start
    assert ct.pos_range.end == t3.pos_range.end
    assert ct.annotated_seq.ref_seq == 'AATAAACCCCCCGGGGGG'
    assert ct.variant_count == len(pam_variants)
    assert ct.annotated_seq.variants == pam_variants
