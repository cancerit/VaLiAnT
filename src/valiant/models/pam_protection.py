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
from itertools import chain
import logging
from typing import Dict, Optional, Set, FrozenSet
import pandas as pd
from pyranges import PyRanges
from pysam import VariantRecord
from .base import GenomicRange
from .sequences import ReferenceSequence
from .variant import BaseVariant, get_variant, SubstitutionVariant
from ..loaders.vcf import get_vcf
from ..utils import get_id_column


@dataclass(frozen=True)
class PamProtectedReferenceSequence(ReferenceSequence):
    __slots__ = {'sequence', 'genomic_range', 'pam_protected_sequence'}

    pam_protected_sequence: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.pam_protected_sequence) != len(self.genomic_range):
            raise ValueError("PAM protected sequence and genomic range have different lengths!")

    @classmethod
    def from_reference_sequence(cls, ref_seq: ReferenceSequence, pam_seq: str) -> PamProtectedReferenceSequence:
        return cls(ref_seq.sequence, ref_seq.genomic_range, pam_seq)

    def get_subsequence(self, genomic_range: GenomicRange) -> PamProtectedReferenceSequence:
        start, end = self.genomic_range.get_relative_subrange(genomic_range)
        ref_seq: ReferenceSequence = super().get_subsequence(genomic_range)
        pam_seq: str = self.pam_protected_sequence[start:end]
        return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_seq)

    def apply_variant(self, variant: BaseVariant, ref_check: bool = False) -> str:
        if not self.genomic_range.contains_position(variant.genomic_position):
            raise ValueError("Variant not in genomic range!")
        # offset: int = variant.genomic_position.position - self.genomic_range.start
        return variant.mutate(self.pam_protected_sequence, self.genomic_range.start, ref_check=ref_check)


@dataclass(frozen=True)
class PamVariant(SubstitutionVariant):
    __slots__ = {'genomic_position', 'ref', 'alt'}

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.ref) != 1 or len(self.alt) != 1:
            raise ValueError(
                "Only single-nucleotide substitutions are allowed "
                "for the purposes of PAM protection!")

    @classmethod
    def from_substitution(cls, var: SubstitutionVariant) -> PamVariant:
        return PamVariant(var.genomic_position, var.ref, var.alt)

    @classmethod
    def from_variant_record(cls, r: VariantRecord) -> PamVariant:
        var: BaseVariant = get_variant(r)
        if not isinstance(var, SubstitutionVariant):
            raise ValueError("PAM protection variants must be substitutions!")
        return cls.from_substitution(var)


class PamProtectionVariantRepository:
    __slots__ = {'sgrna_ids', '_variants', '_ranges'}

    def __init__(self, sgrna_ids: Set[str] = None) -> None:
        self.sgrna_ids = sgrna_ids or set()
        self._variants: Dict[str, Set[PamVariant]] = {
            sgrna_id: set()
            for sgrna_id in sgrna_ids
        } if sgrna_ids else {}
        self._ranges: Optional[PyRanges] = None

    @property
    def count(self) -> int:
        return sum(map(len, self._variants.values()))

    def get_sgrna_variants(self, sgrna_id: str) -> FrozenSet[PamVariant]:
        return frozenset(self._variants[sgrna_id])

    def get_sgrna_variants_bulk(self, sgrna_ids: FrozenSet[str]) -> FrozenSet[PamVariant]:
        try:
            return frozenset(chain.from_iterable(
                self.get_sgrna_variants(sgrna_id)
                for sgrna_id in sgrna_ids
            ))

        except KeyError as ex:
            sgrna_id: str = ex.args[0]
            raise RuntimeError(f"sgRNA ID '{sgrna_id} not found!'")

    def load(self, fp: str) -> None:

        # Load variants from VCF
        with get_vcf(fp) as variant_file:
            for record in variant_file.fetch():
                sgrna_id: str = record.info['SGRNA'].strip()
                if sgrna_id in self._variants:
                    self._variants[sgrna_id].add(
                        PamVariant.from_variant_record(record))

        # Log loaded variant statistics
        logging.debug("Collected %d PAM protection variants." % self.count)
        for sgrna_id in self.sgrna_ids:
            if not self._variants[sgrna_id]:
                logging.info("No PAM protection variants for sgRNA %s." % sgrna_id)

        # Populate genomic range table
        df: pd.DataFrame = pd.DataFrame.from_records([
            (*variant.get_pyrange_record(), sgrna_id)
            for sgrna_id, variants in self._variants.items()
            for variant in variants
        ], columns=[
            'Chromosome',
            'Start',
            'End',
            'sgrna_id'
        ])
        df.sgrna_id = df.sgrna_id.astype('category')
        df['variant_id'] = get_id_column(df.shape[0])
        self._ranges = PyRanges(df)


def compute_pam_protected_sequence(
    ref_seq: ReferenceSequence,
    pam_variants: Set[PamVariant]
) -> PamProtectedReferenceSequence:
    pam_seq: str = ref_seq.sequence

    if len(pam_variants) > 0:
        for variant in pam_variants:

            # Validate variant genomic position relative to the sequence's
            offset: int = variant.get_ref_offset(ref_seq)

            # Update PAM-protected sequence
            pam_seq = variant.mutate_from(pam_seq, offset, ref_check=True)

    return PamProtectedReferenceSequence.from_reference_sequence(ref_seq, pam_seq)
