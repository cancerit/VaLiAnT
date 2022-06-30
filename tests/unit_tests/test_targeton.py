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
    genomic_range =  GenomicRange(chromosome, 500, 500 + len(ref_seq) - 1, '+')
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

        # Verify the temporay field was dropped

        df = mutation_collections[mutator].df
        assert '_codon_index' not in df.columns

        # Validate the sgRNA ID field
        df['codon_index'] = (df['mut_position'] + frame) // 3
        mask = df['codon_index'].mod(2).ne(0)
        assert pd.isna(df.loc[mask, 'meta_pam_mut_sgrna_id'].unique())
        assert df.loc[df.codon_index == 0, 'meta_pam_mut_sgrna_id'].unique() == sgrna_id_1
        assert df.loc[df.codon_index == 2, 'meta_pam_mut_sgrna_id'].unique() == sgrna_id_2
