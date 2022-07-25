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

import pandas as pd
import pytest
from valiant.enums import VariantType
from valiant.models.base import GenomicPosition
from valiant.models.custom_variants import CustomVariantMutation, CustomVariantMutationCollection
from valiant.models.variant import CustomVariant, DeletionVariant, InsertionVariant, SubstitutionVariant

deletion_variant = DeletionVariant(GenomicPosition('X', 5), 'A')
insertion_variant = InsertionVariant(GenomicPosition('X', 5), 'C')
substitution_variant = SubstitutionVariant(GenomicPosition('X', 5), 'A', 'C')

variants = [
    deletion_variant,
    insertion_variant,
    substitution_variant
]

variant_mseqs = [
    (deletion_variant, 'AAAAAAAAA'),
    (insertion_variant, 'AAAACAAAAAA'),
    (substitution_variant, 'AAAACAAAAA')
]


@pytest.mark.parametrize('variant', variants)
def test_custom_variant_mutation_init(variant):
    seq = 'AAAAAAAAAA'
    custom_variant = CustomVariant(variant, 'vcf_alias', 'VARIANT_ID')

    # Initialise custom variant mutation
    CustomVariantMutation(custom_variant, seq, False, frozenset())


@pytest.mark.parametrize('variant', variants)
def test_custom_variant_mutation_to_row(variant):
    seq = 'AAAAAAAAAA'
    var_type = variant.type
    custom_variant = CustomVariant(variant, 'vcf_alias', 'VARIANT_ID')

    # Initialise custom variant mutation
    cvm = CustomVariantMutation(custom_variant, seq, False, frozenset())

    # Check custom variant mutation
    vcf_alias, vcf_variant_id, var_type_, pos, ref, alt, seq_, sgrna_ids, in_const = cvm.to_row()
    assert vcf_alias == custom_variant.vcf_alias
    assert vcf_variant_id == custom_variant.vcf_variant_id
    assert var_type_ == var_type.value
    assert pos == custom_variant.base_variant.genomic_position.position
    assert sgrna_ids == ''
    assert in_const == 0

    if var_type != VariantType.INSERTION:
        assert ref == custom_variant.base_variant.ref
    else:
        assert ref is None

    if var_type != VariantType.DELETION:
        assert alt == custom_variant.base_variant.alt
    else:
        assert alt is None

    assert seq_ == seq


@pytest.mark.parametrize('variant,mseq', variant_mseqs)
def test_custom_variant_mutation_collection_init(variant, mseq):
    var_type = variant.type
    var_pos = variant.genomic_position.position
    var_ref = getattr(variant, 'ref', None)
    var_alt = getattr(variant, 'alt', None)
    df = pd.DataFrame.from_records([
        ('VCF', 'ID', var_type.value, var_pos, var_ref, var_alt, mseq)
    ])

    # Initialise mutation collection
    cvmc = CustomVariantMutationCollection(df)

    # Check data frame
    assert cvmc.df.equals(df)


@pytest.mark.parametrize('variant,mseq', variant_mseqs)
def test_custom_variant_mutation_collection_from_variants(variant, mseq):
    vcf_alias = 'vcf_alias'
    vcf_var_id = 'VARIANT_ID'
    var_type = variant.type
    custom_variant = CustomVariant(variant, vcf_alias, vcf_var_id)
    cvm = CustomVariantMutation(custom_variant, mseq, False, frozenset())

    # Initialise mutation collection
    cvmc = CustomVariantMutationCollection.from_variants([cvm])

    # Check data frame
    assert cvmc.df.shape[0] == 1
    var_row = cvmc.df.iloc[0]
    assert var_row.vcf_alias == vcf_alias
    assert var_row.vcf_var_id == vcf_var_id
    assert var_row.mut_position == variant.genomic_position.position
    assert var_row.mseq == mseq

    if var_type != VariantType.INSERTION:
        assert var_row.ref == variant.ref
    else:
        assert pd.isna(var_row.ref)

    if var_type != VariantType.DELETION:
        assert var_row.new == variant.alt
    else:
        assert pd.isna(var_row.new)
