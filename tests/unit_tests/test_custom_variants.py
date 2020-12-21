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
    CustomVariantMutation(custom_variant, seq)


@pytest.mark.parametrize('variant', variants)
def test_custom_variant_mutation_to_row(variant):
    seq = 'AAAAAAAAAA'
    var_type = variant.type
    custom_variant = CustomVariant(variant, 'vcf_alias', 'VARIANT_ID')

    # Initialise custom variant mutation
    cvm = CustomVariantMutation(custom_variant, seq)

    # Check custom variant mutation
    vcf_alias, vcf_variant_id, var_type_, pos, ref, alt, seq_ = cvm.to_row()
    assert vcf_alias == custom_variant.vcf_alias
    assert vcf_variant_id == custom_variant.vcf_variant_id
    assert var_type_ == var_type.value
    assert pos == custom_variant.base_variant.genomic_position.position

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
    cvm = CustomVariantMutation(custom_variant, mseq)

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
