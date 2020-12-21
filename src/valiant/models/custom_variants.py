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

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd
from .mutated_sequences import BaseMutationCollection
from .oligo_renderer import BaseOligoRenderer
from .pam_protection import PamProtectedReferenceSequence
from .variant import BaseVariant, CustomVariant


@dataclass(init=False)
class CustomVariantOligoRenderer(BaseOligoRenderer):
    __slots__ = {
        'ref_seq',
        'gene_id',
        'transcript_id',
        'adaptor_5',
        'adaptor_3',
        '_oligo_name_prefix'
    }

    def __init__(
        self,
        ref_seq: PamProtectedReferenceSequence,
        gene_id: str,
        transcript_id: str,
        adaptor_5: str,
        adaptor_3: str
    ) -> None:
        super().__init__(ref_seq, gene_id, transcript_id, adaptor_5, adaptor_3)

    def get_oligo_name(self, vcf_alias: str, var_type: int, start: int, ref: Optional[str], alt: Optional[str]) -> str:
        return super()._get_oligo_name(var_type, vcf_alias, start, ref, alt)


@dataclass
class CustomVariantMutation:
    variant: CustomVariant
    sequence: str

    def to_row(self) -> Tuple[Optional[str], Optional[str], int, int, Optional[str], Optional[str], str]:
        var: BaseVariant = self.variant.base_variant
        return (
            self.variant.vcf_alias,
            self.variant.vcf_variant_id,
            var.type.value,
            var.genomic_position.position,
            getattr(var, 'ref', None),
            getattr(var, 'alt', None),
            self.sequence
        )


@dataclass
class CustomVariantMutationCollection(BaseMutationCollection):
    __slots__ = {'df'}

    @classmethod
    def from_variants(cls, mutations: List[CustomVariantMutation]) -> CustomVariantMutationCollection:
        df: pd.DataFrame = pd.DataFrame.from_records([
            mutation.to_row()
            for mutation in mutations
        ], columns=[
            'vcf_alias',
            'vcf_var_id',
            'var_type',
            'mut_position',
            'ref',
            'new',
            'mseq'
        ])

        # Compress table
        df.vcf_alias = df.vcf_alias.astype('category')
        df.vcf_var_id = df.vcf_var_id.astype('string')
        df.var_type = df.var_type.astype(np.int8)
        df.mut_position = df.mut_position.astype(np.int32)
        df.ref = df.ref.astype('category')
        df.new = df.new.astype('category')
        df.mseq = df.mseq.astype('string')

        return cls(df=df)
