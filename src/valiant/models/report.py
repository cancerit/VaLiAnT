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

from __future__ import annotations
from dataclasses import dataclass
import json
from typing import Dict, Any
from .. import __version__ as APP_VERSION, __package__ as PACKAGE_NAME
from .options import Options


@dataclass
class Report:
    min_length: int
    max_length: int
    rc: bool

    @classmethod
    def from_options(cls, options: Options) -> Report:
        return cls(
            min_length=options.oligo_min_length,
            max_length=options.oligo_max_length,
            rc=options.revcomp_minus_strand
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            'appName': PACKAGE_NAME,
            'appVersion': APP_VERSION,
            'reverseComplementOnMinusStrand': self.rc,
            'minOligoLength': self.min_length,
            'maxOligoLength': self.max_length
        }

    def write(self, fp: str) -> None:
        with open(fp, 'w') as fh:
            json.dump(self.to_dict(), fh, separators=(',', ':'))
