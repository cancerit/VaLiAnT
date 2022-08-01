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

import abc
import logging
from typing import Optional

from ..errors import InvalidConfig
from ..utils import is_adaptor_valid
from .options import Options
from pydantic import BaseModel


class BaseConfig(BaseModel, abc.ABC):

    # Metadata
    species: str
    assembly: str

    # Adaptors
    adaptor_5: Optional[str]
    adaptor_3: Optional[str]

    # Filters
    min_length: int
    max_length: int

    # Paths
    codon_table_fp: Optional[str]
    oligo_info_fp: str
    ref_fasta_fp: str
    output_dir: str

    class Config:
        allow_population_by_field_name = True
        fields = {
            'species': 'species',
            'assembly': 'assembly',
            'adaptor_5': 'adaptor5',
            'adaptor_3': 'adaptor3',
            'min_length': 'minOligoLength',
            'max_length': 'maxOligoLength',
            'codon_table_fp': 'codonTableFilePath',
            'oligo_info_fp': 'oligoInfoFilePath',
            'ref_fasta_fp': 'refFASTAFilePath',
            'output_dir': 'outputDirPath'
        }

    def __post_init__(self) -> None:
        if not self.is_valid():
            raise InvalidConfig()

    @abc.abstractmethod
    def get_options(self) -> Options:
        pass

    def write(self, fp: str) -> None:
        with open(fp, 'w') as fh:
            fh.write(self.json(by_alias=True, separators=(',', ':')))

    def is_valid(self) -> bool:
        success: bool = True

        # Validate adaptors
        for adaptor in [self.adaptor_5, self.adaptor_3]:
            if not is_adaptor_valid(adaptor):
                logging.critical("Invalid adaptor sequence '%s'!" % adaptor)
                success = False

        return success
