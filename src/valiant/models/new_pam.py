########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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
from typing import List

from .new_seq import AltSeqBuilder
from .pam_protection import PamVariant
from .variant import BaseVariantT


@dataclass(frozen=True, init=False)
class PamBgAltSeqBuilder(AltSeqBuilder):
    def __init__(
        self,
        start: int,
        sequence: str,
        bg_variants: List[BaseVariantT],
        pam_variants: List[PamVariant],
        cds_prefix: str = '',
        cds_suffix: str = ''
    ) -> None:
        if len(cds_suffix) > 2 or len(cds_prefix) > 2:
            raise ValueError("Invalid CDS extension!")

        super().__init__(
            start,
            sequence, [
                bg_variants,
                pam_variants
            ],
            prefix=cds_prefix,
            suffix=cds_suffix)
