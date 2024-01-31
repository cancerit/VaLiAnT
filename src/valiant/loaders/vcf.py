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
from contextlib import contextmanager, closing
from functools import partial
from typing import Callable, Generator, TypeVar

from pysam import VariantRecord, VariantFile

from ..contig_filter import ContigFilter
from ..variant import VariantWithContigT


VcfId = str | None
GetVcfId = Callable[[VariantRecord], VcfId]
GetVcfVar = Callable[[VcfId, VariantRecord], VariantWithContigT]
FilterVar = Callable[[VariantWithContigT], bool]


@contextmanager
def open_vcf(fp: str, **kwargs) -> Generator[VariantFile, None, None]:
    with closing(VariantFile(fp, **kwargs)) as vcf:
        yield vcf


def get_var_id_from_info(vcf_id_tag: str, r: VariantRecord) -> str | None:
    vcf_var_id = r.info.get(vcf_id_tag, None)
    return str(vcf_var_id) if vcf_var_id is not None else None


def get_var_id_from_id(r: VariantRecord) -> str | None:
    return str(r.id) if r.id is not None else None


def load_vcf(
    fp: str,
    ft: ContigFilter,
    parse_f: GetVcfVar,
    filter_f: FilterVar | None = None,
    vcf_id_tag: str | None = None
) -> list[VariantWithContigT]:

    if not vcf_id_tag:
        logging.info("No INFO tag specified as variant identifier for VCF file '%s', falling back to the ID field." % fp)

    # Set the function to extract the record identifier
    get_id_f = partial(get_var_id_from_info, vcf_id_tag) if vcf_id_tag else get_var_id_from_id

    def vcf_record_to_variant(r: VariantRecord) -> VariantWithContigT:
        return parse_f(get_id_f(r), r)

    def filter_record(r: VariantRecord) -> bool:
        # Evaluate the contig
        return r.contig in ft.contigs

    def filter_variant(v: VariantWithContigT) -> bool:
        # Evaluate the corrected position for insertions and deletions
        # TODO: consider if it's worth evaluating the end position as well
        if v.pos not in ft.contigs[v.contig]:
            return False
        return filter_f(v) if filter_f is not None else True

    # TODO: if indexed, fetch by chromosome?
    #  Skip record filtering in that case.
    with open_vcf(fp) as vcf:
        if vcf_id_tag and vcf_id_tag not in vcf.header.info:
            raise ValueError(f"INFO tag '{vcf_id_tag}' not found in VCF file '{fp}'!")

        return [
            t for t in (
                vcf_record_to_variant(r)
                for r in (
                    record
                    for record in vcf.fetch()
                    if filter_record(record)
                )
            )
            if filter_variant(t)
        ]
