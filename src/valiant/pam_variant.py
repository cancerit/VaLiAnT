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

from __future__ import annotations

from dataclasses import dataclass

from pysam import VariantRecord

from .contig_filter import ContigFilter
from .loaders.vcf import load_vcf
from .variant import VariantWithContig

from .custom_variant import CustomVariant


VCF_SGRNA_ID_INFO = 'SGRNA'


class InvalidPamVariant(Exception):
    pass


@dataclass
class PamVariant(VariantWithContig):
    sgrna_id: str

    def __post_init__(self) -> None:
        assert self.ref_len == 1 and self.alt_len == 1

    @classmethod
    def from_record_with_id(cls, id: str | None, r: VariantRecord) -> PamVariant:
        var = CustomVariant.from_record_with_id(id, r)
        try:
            sgrna_id = r.info[VCF_SGRNA_ID_INFO].strip()
        except KeyError:
            raise InvalidPamVariant("PAM protection edit with no sgRNA ID!")
        return cls(var.pos, var.ref, var.alt, var.contig, sgrna_id)

    @classmethod
    def load_vcf(cls, fp: str, ft: ContigFilter, vcf_id_tag: str | None = None) -> list[PamVariant]:
        return load_vcf(fp, ft, cls.from_record_with_id, vcf_id_tag=vcf_id_tag)
