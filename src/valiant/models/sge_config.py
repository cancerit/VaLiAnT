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

from dataclasses import dataclass
from typing import Any, ClassVar, Dict, Optional
from ..enums import ExecMode
from .config import BaseConfig
from .options import Options


@dataclass
class SGEConfig(BaseConfig):
    mode: ClassVar[ExecMode] = ExecMode.SGE

    revcomp_minus_strand: bool
    gff_fp: Optional[str]
    pam_fp: Optional[str]
    vcf_fp: Optional[str]

    def get_options(self) -> Options:
        return Options(
            revcomp_minus_strand=self.revcomp_minus_strand,
            oligo_max_length=self.max_length,
            oligo_min_length=self.min_length)

    def to_dict(self) -> Dict[str, Any]:
        return {
            **super().to_dict(),
            'GFFFilePath': self.gff_fp,
            'PAMProtectionVCFFilePath': self.pam_fp,
            'customVCFManifestFilePath': self.vcf_fp
        }
