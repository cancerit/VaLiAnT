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

from __future__ import annotations

import logging
from dataclasses import dataclass

from .csv import load_csv
from ..strings.file_path import FilePath
from ..utils import has_duplicates, get_dataclass_fields


@dataclass(slots=True)
class VcfInfo:
    # BEWARE: field names and order must be preserved
    #  as they are used to validate the header and parse the CSV rows
    vcf_alias: str
    vcf_id_tag: str | None
    vcf_path: FilePath

    def __post_init__(self) -> None:
        assert self.vcf_alias
        assert self.vcf_path

    @classmethod
    def from_list(cls, a: list[str | None]) -> VcfInfo:
        return cls(a[0], a[1], FilePath(a[2]))


@dataclass(slots=True)
class VcfManifest:
    vcfs: list[VcfInfo]

    def __post_init__(self) -> None:
        aliases = [x.vcf_alias for x in self.vcfs]
        if has_duplicates(aliases):
            raise ValueError("Duplicate VCF file aliases found in manifest!")
        paths = [x.vcf_path for x in self.vcfs]
        if has_duplicates(paths):
            raise ValueError("Duplicate VCF file paths found in manifest!")

    @classmethod
    def load(cls, fp: str) -> VcfManifest:
        # TODO: catch header validation failure
        # TODO: consider whether to allow any column order as in v. 3.*
        vcf_infos = [
            VcfInfo.from_list(ls)
            for ls in load_csv(fp, get_dataclass_fields(VcfInfo))
        ]
        logging.debug(f"{len(vcf_infos)} VCF files listed in manifest.")
        return cls(vcf_infos)

    def test_vcf_extance(self) -> None:
        for vcf_info in self.vcfs:
            vcf_info.vcf_path.test_extance()
