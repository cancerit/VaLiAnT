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

from typing import List, Optional
from .config import BaseConfig
from .options import Options


class CDNAConfig(BaseConfig):
    annot_fp: Optional[str]

    class Config:
        fields = {'annot_fp': 'annotationFilePath'}

    def get_options(self) -> Options:
        return Options(
            revcomp_minus_strand=False,
            oligo_max_length=self.max_length,
            oligo_min_length=self.min_length)

    @property
    def input_file_paths(self) -> List[str]:
        fps = super().input_file_paths
        if self.annot_fp is not None:
            fps.append(self.annot_fp)
        return fps
