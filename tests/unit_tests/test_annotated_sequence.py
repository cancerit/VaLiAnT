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
import pytest
from valiant.models.annotated_sequence import AnnotatedSequencePair, CDSAnnotatedSequencePair
from valiant.models.base import GenomicPosition, StrandedPositionRange
from valiant.models.variant import SubstitutionVariant
from .utils import load_codon_table


codon_table = load_codon_table()


def test_get_variant_mutation_types():

    # Generate annotated sequence pair
    genomic_positions = [100, 101]
    sp = AnnotatedSequencePair(
        StrandedPositionRange(100, 105, '+'),
        'AAGTACG',
        'CCGTACG', [
            SubstitutionVariant(GenomicPosition('X', pos), 'A', 'C')
            for pos in genomic_positions
        ])

    print(sp)

    # Retrieve mutation types
    mutation_types = sp.get_variant_mutation_types(codon_table)

    assert len(mutation_types) == len(genomic_positions)


@pytest.mark.parametrize('no_duplicate_codons', [True, False])
def test_cds_get_variant_mutation_types(no_duplicate_codons):

    # Generate annotated sequence pair
    genomic_positions = [100, 101]
    sp = CDSAnnotatedSequencePair(
        StrandedPositionRange(100, 105, '+'),
        'AAGTACG',
        'CCGTACG', [
            SubstitutionVariant(GenomicPosition('X', pos), 'A', 'C')
            for pos in genomic_positions
        ],
        'T',
        'C')

    # Retrieve mutation types
    with pytest.raises(ValueError) if no_duplicate_codons else nullcontext():
        mutation_types = sp.get_variant_mutation_types(
            codon_table, no_duplicate_codons=no_duplicate_codons)

        assert len(mutation_types) == len(genomic_positions)
