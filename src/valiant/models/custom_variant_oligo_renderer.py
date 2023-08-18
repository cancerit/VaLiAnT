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

from dataclasses import dataclass
from typing import Optional

from .new_pam import PamBgAltSeqBuilder
from .oligo_renderer import BaseOligoRenderer


@dataclass(init=False)
class CustomVariantOligoRenderer(BaseOligoRenderer):
    __slots__ = {
        'ref_seq',
        'gene_id',
        'transcript_id',
        'adaptor_5',
        'adaptor_3',
        '_oligo_name_prefix'
    }

    def __init__(
        self,
        ref_seq: PamBgAltSeqBuilder,
        gene_id: str,
        transcript_id: str,
        adaptor_5: str,
        adaptor_3: str
    ) -> None:
        super().__init__(ref_seq, gene_id, transcript_id, adaptor_5, adaptor_3)

    def get_oligo_name(self, vcf_alias: str, var_type: int, start: int, ref: Optional[str], alt: Optional[str]) -> str:
        return super()._get_oligo_name(var_type, vcf_alias, start, ref, alt)
