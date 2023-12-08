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

import abc
import logging
import os
from typing import Any, List, Optional

from pydantic import BaseModel, Field

from .errors import InvalidConfig
from .options import Options
from .utils import is_adaptor_valid


class BaseConfig(BaseModel, abc.ABC):

    # Metadata
    species: str = Field()
    assembly: str = Field()

    # Adaptors
    adaptor_5: Optional[str] = Field(alias='adaptor5')
    adaptor_3: Optional[str] = Field(alias='adaptor3')

    # Filters
    min_length: int = Field(alias='minOligoLength')
    max_length: int = Field(alias='maxOligoLength')

    # Paths
    codon_table_fp: str | None = Field(alias='codonTableFilePath')
    oligo_info_fp: str = Field(alias='oligoInfoFilePath')
    ref_fasta_fp: str = Field(alias='refFASTAFilePath')
    output_dir: str = Field(alias='outputDirPath')

    class Config:
        populate_by_name = True

    def __init__(__pydantic_self__, **data: Any) -> None:
        super().__init__(**data)
        if not __pydantic_self__.is_valid():
            raise InvalidConfig()

    @abc.abstractmethod
    def get_options(self) -> Options:
        pass

    @property
    def input_file_paths(self) -> List[str]:
        fps: List[str] = [
            self.oligo_info_fp,
            self.ref_fasta_fp
        ]
        if self.codon_table_fp is not None:
            fps.append(self.codon_table_fp)
        return fps

    def write(self, fp: str) -> None:
        with open(fp, 'w') as fh:
            fh.write(self.model_dump_json(by_alias=True))

    def is_valid(self) -> bool:
        success: bool = True

        # Validate adaptors
        for adaptor in [self.adaptor_5, self.adaptor_3]:
            if not is_adaptor_valid(adaptor):
                logging.error("Invalid adaptor sequence '%s'!" % adaptor)
                success = False

        # Validate oligonucleotide lengths
        for length, label in [
            (self.max_length, "maximum"),
            (self.min_length, "minimum"),
        ]:
            if length < 1:
                logging.error(
                    "Invalid %s oligonucleotide length: "
                    "not strictly positive!" % label)
                success = False

        return success

    def get_output_file_path(self, fp: str) -> str:
        return os.path.join(self.output_dir, fp)
