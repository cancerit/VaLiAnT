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
import os
import pandas as pd
from ..constants import METADATA_FIELDS_SET, METADATA_FIELDS
from ..loaders.vcf import write_vcf
from ..utils import get_constant_category
from ..writers import write_oligo_metadata, write_oligo_unique
from .oligo_generation_info import OligoGenerationInfo
from .refseq_repository import ReferenceSequenceRepository
from .variant import get_records


def _fill_metadata(species: str, assembly: str, metadata: pd.DataFrame) -> pd.DataFrame:

    # Add global metadata
    rown: int = metadata.shape[0]
    metadata['species'] = get_constant_category(species, rown)
    metadata['assembly'] = get_constant_category(assembly, rown)

    # Add missing columns
    for field in METADATA_FIELDS_SET - set(metadata.columns):
        metadata[field] = None

    return metadata


@dataclass(frozen=True)
class MetadataTable:
    __slots__ = {'metadata', 'oligo_length_mask', 'short_oligo_n', 'long_oligo_n'}

    metadata: pd.DataFrame
    oligo_length_mask: pd.Series
    short_oligo_n: int
    long_oligo_n: int

    def __init__(self, metadata: pd.DataFrame, oligo_max_length: int) -> None:
        if oligo_max_length < 1:
            raise ValueError("Invalid maximum oligonucleotide length!")

        object.__setattr__(self, 'metadata', metadata)
        object.__setattr__(self, 'oligo_length_mask', metadata.oligo_length <= oligo_max_length)
        object.__setattr__(self, 'short_oligo_n', self.oligo_length_mask.sum())
        object.__setattr__(self, 'long_oligo_n', len(self.oligo_length_mask) - self.short_oligo_n)

    @classmethod
    def from_partial(
        cls,
        species: str,
        assembly: str,
        metadata: pd.DataFrame,
        oligo_max_length: int
    ) -> MetadataTable:
        return cls(
            _fill_metadata(species, assembly, metadata),
            oligo_max_length)

    def write_metadata_file(self, fp: str) -> None:
        write_oligo_metadata((
            self.metadata.loc[self.oligo_length_mask, METADATA_FIELDS] if self.long_oligo_n > 0 else
            self.metadata[METADATA_FIELDS]
        ), fp)

    def write_excluded_metadata_file(self, fp: str) -> None:
        write_oligo_metadata(
            self.metadata.loc[~self.oligo_length_mask, METADATA_FIELDS], fp)

    def write_unique_file(self, fp: str) -> None:
        write_oligo_unique(
            self.metadata.loc[
                self.oligo_length_mask,
                ['oligo_name', 'mseq']
            ].sort_values(
                'oligo_name',
                kind='mergesort'
            ).drop_duplicates(subset=['mseq'], keep='first'),
            fp)

    def write_vcf_file(self, fp: str, ref_repository: ReferenceSequenceRepository) -> None:
        metadata = (
            self.metadata[self.oligo_length_mask] if self.long_oligo_n > 0 else
            self.metadata
        )
        write_vcf(
            fp,
            list(metadata.ref_chr.cat.categories.values),
            get_records(ref_repository, metadata))

    def write_common_files(self, out_dir: str, base_fn: str) -> None:
        if self.short_oligo_n > 0:

            # Save metadata to file (filtered and reordered columns)
            metadata_fn = base_fn + '_meta.csv'
            self.write_metadata_file(os.path.join(out_dir, metadata_fn))

            # Save unique oligonucleotides to file (and filter by length)
            unique_oligos_fn = base_fn + '_unique.csv'
            self.write_unique_file(os.path.join(out_dir, unique_oligos_fn))

        if self.long_oligo_n > 0:

            # Save discarded metadata to file
            excluded_metadata_fn = base_fn + '_meta_excluded.csv'
            self.write_excluded_metadata_file(
                os.path.join(out_dir, excluded_metadata_fn))

    def write_sge_files(
        self,
        out_dir: str,
        base_fn: str,
        ref_repository: ReferenceSequenceRepository
    ) -> None:
        self.write_common_files(out_dir, base_fn)
        if self.short_oligo_n > 0:

            # Save variants to file (VCF format)
            vcf_fn = base_fn + '.vcf'
            self.write_vcf_file(os.path.join(out_dir, vcf_fn), ref_repository)

    def get_info(self) -> OligoGenerationInfo:
        return OligoGenerationInfo(self.long_oligo_n)
