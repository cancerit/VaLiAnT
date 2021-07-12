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
import pytest
from valiant.models.variant import SubstitutionVariant, DeletionVariant, InsertionVariant, get_variant_from_tuple
from valiant.models.base import GenomicPosition


POSITION = GenomicPosition('X', 10000)


@pytest.mark.parametrize('pos,ref,alt,variant_type', [
    (3, 'C', 'CGTTAAGGTGTCAAGGAAAGGAA', InsertionVariant),
    (1, 'A', 'GTTAAGGTGTCAAGGAAAGGAAA', InsertionVariant),
    (1, 'A', 'GTTAAGGTGTCAAGGAAAGGAAC', SubstitutionVariant),
    (3, 'CGTTAAGGTGTCAAGGAAAGGAA', 'C', DeletionVariant),
    (1, 'GTTAAGGTGTCAAGGAAAGGAAA', 'A', DeletionVariant),
    (3, 'A', 'C', SubstitutionVariant),
    (1, 'G', 'T', SubstitutionVariant)
])
def test_get_variant_from_tuple(pos, ref, alt, variant_type):
    valid = variant_type is not None
    with pytest.raises(ValueError) if not valid else nullcontext():
        variant = get_variant_from_tuple('X', pos, ref, alt)
        assert isinstance(variant, variant_type)


@pytest.mark.parametrize('seq,offset,ref,alt,mseq', [
    ('AAACCCGGGTTT', 3, 'CCC', 'TTT', 'AAATTTGGGTTT')
])
def test_substitution_variant_mutate(seq, offset, ref, alt, mseq):
    var = SubstitutionVariant(POSITION, ref, alt)
    assert var.mutate_from(seq, offset) == mseq


@pytest.mark.parametrize('seq,offset,ref,mseq', [
    ('AAACCCGGGTTT', 3, 'CCC', 'AAAGGGTTT')
])
def test_deletion_variant_mutate(seq, offset, ref, mseq):
    var = DeletionVariant(POSITION, ref)
    assert var.mutate_from(seq, offset) == mseq


@pytest.mark.parametrize('seq,offset,alt,mseq', [
    ('AAACCCGGGTTT', 3, 'TTT', 'AAATTTCCCGGGTTT')
])
def test_insertion_variant_mutate(seq, offset, alt, mseq):
    var = InsertionVariant(POSITION, alt)
    assert var.mutate_from(seq, offset) == mseq
