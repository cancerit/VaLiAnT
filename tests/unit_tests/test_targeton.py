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

import pytest
from valiant.enums import MutationType
from valiant.models.annotated_sequence import CDSAnnotatedSequencePair
from valiant.models.base import GenomicPosition, GenomicRange
from valiant.models.pam_protection import PamVariant
from valiant.models.targeton import PamProtCDSTargeton
from .utils import load_codon_table


SGRNA_ID = 'some-id'
codon_table = load_codon_table()


def test_pam_prot_targeton_from_reference_sequence():
    chromosome = 'X'
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

    chromosome = 'X'
    PamProtCDSTargeton(CDSAnnotatedSequencePair(
        GenomicRange(chromosome, 500, 505, '+'),
        'GCGTGC',
        f"G{alt}GTGC",
        [PamVariant(GenomicPosition(chromosome, 501), 'C', alt, SGRNA_ID)],
        'T',
        'CA'))
