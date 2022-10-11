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


def get_default_cds_annotate_sequence_pair():
    spr = StrandedPositionRange(100, 106, '+')
    genomic_positions = [100, 101]
    return CDSAnnotatedSequencePair(
        spr,
        'AAGTACG',
        'CCGTACG', [
            SubstitutionVariant(GenomicPosition('X', pos), 'A', 'C')
            for pos in genomic_positions
        ],
        'T',
        'C')


cds_sequence_pair = get_default_cds_annotate_sequence_pair()
last_codon_index = (cds_sequence_pair.ext_seq_length // 3) - 1


def test_get_variant_mutation_types():

    # Retrieve mutation types
    mutation_types = cds_sequence_pair.get_variant_mutation_types(codon_table)

    assert len(mutation_types) == cds_sequence_pair.variant_count


@pytest.mark.parametrize('no_duplicate_codons', [True, False])
def test_cds_get_variant_mutation_types(no_duplicate_codons):

    # Retrieve mutation types
    with pytest.raises(ValueError) if no_duplicate_codons else nullcontext():
        mutation_types = cds_sequence_pair.get_variant_mutation_types(
            codon_table, no_duplicate_codons=no_duplicate_codons)

        assert len(mutation_types) == cds_sequence_pair.variant_count


# Position ranges assuming the reading frame is one
@pytest.mark.parametrize('spr,codon_indices', [
    (StrandedPositionRange(start, end, '+'), codon_indices)
    for start, end, codon_indices in [
        (
            cds_sequence_pair.pos_range.start,
            cds_sequence_pair.pos_range.end,
            list(range(cds_sequence_pair.ext_seq_length // 3))
        ),
        (
            cds_sequence_pair.pos_range.start,
            cds_sequence_pair.pos_range.start,
            [0]
        ),
        (
            cds_sequence_pair.pos_range.start,
            cds_sequence_pair.pos_range.start + 1,
            [0]
        ),
        (
            cds_sequence_pair.pos_range.start,
            cds_sequence_pair.pos_range.start + 2,
            [0, 1]
        ),
        (
            cds_sequence_pair.pos_range.end,
            cds_sequence_pair.pos_range.end,
            [last_codon_index]
        ),
        (
            cds_sequence_pair.pos_range.end - 1,
            cds_sequence_pair.pos_range.end,
            [last_codon_index]
        ),
        (
            cds_sequence_pair.pos_range.end - 2,
            cds_sequence_pair.pos_range.end,
            [last_codon_index - 1, last_codon_index]
        )
    ]
])
def test_get_codon_indices_in_range(spr, codon_indices):

    # Generate annotated sequence pair
    assert cds_sequence_pair.get_codon_indices_in_range(spr) == codon_indices
