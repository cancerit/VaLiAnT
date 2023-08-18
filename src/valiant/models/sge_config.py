########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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

import logging
from typing import List, Optional

from .config import BaseConfig
from .options import Options
from ..utils import get_not_none


class SGEConfig(BaseConfig):
    revcomp_minus_strand: bool

    # File paths
    gff_fp: Optional[str]
    bg_fp: Optional[str]
    pam_fp: Optional[str]
    vcf_fp: Optional[str]

    # Background variants
    force_bg_ns: bool
    force_bg_fs: bool

    class Config:
        fields = {
            'revcomp_minus_strand': 'reverseComplementOnMinusStrand',
            'gff_fp': 'GFFFilePath',
            'bg_fp': 'backgroundVCFFilePath',
            'pam_fp': 'PAMProtectionVCFFilePath',
            'vcf_fp': 'customVCFManifestFilePath',
            'force_bg_ns': 'forceBackgroundNonSynonymous',
            'force_bg_fs': 'forceBackgroundFrameShifting'
        }

    @property
    def input_file_paths(self) -> List[str]:
        return super().input_file_paths + get_not_none([
            self.gff_fp,
            self.bg_fp,
            self.pam_fp,
            self.vcf_fp
        ])

    def get_options(self) -> Options:
        return Options(
            revcomp_minus_strand=self.revcomp_minus_strand,
            oligo_max_length=self.max_length,
            oligo_min_length=self.min_length)

    def is_valid(self) -> bool:
        success = super().is_valid()
        if self.force_bg_fs and not self.force_bg_ns:
            logging.error(
                "For frame-shifting background variants to be allowed, "
                "non-synonymous variants must also be allowed!")
            success = False
        return success
