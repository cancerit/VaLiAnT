########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023, 2024 Genome Research Ltd
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

from contextlib import contextmanager
from typing import Generator

from pysam import VariantFile, VariantHeader

from .loaders.vcf import open_vcf


def _get_vcf_header_info_items(id_: str, type_: str, n: int) -> list[tuple[str, str | int]]:
    return [
        ('ID', id_),
        ('Type', type_),
        ('Number', n)
    ]


VCF_HEADER_INFO_ITEMS: list[list[tuple[str, str | int]]] = [
    _get_vcf_header_info_items(id_, type_, n)
    for id_, type_, n in [
        ('SGE_SRC', 'String', 1),
        ('SGE_REF', 'String', 1),
        ('SGE_OLIGO', 'String', 1),
        ('SGE_VCF_ALIAS', 'String', 1),
        ('SGE_VCF_VAR_ID', 'String', 1)
    ]
]


@contextmanager
def _open_vcf_write(fp, header: VariantHeader) -> Generator[VariantFile, None, None]:
    with open_vcf(fp, mode='w', header=header) as vcf:
        yield vcf


def _get_vcf_header(contigs: list[str]) -> VariantHeader:
    header: VariantHeader = VariantHeader()
    header.add_meta('source', value='valiant')

    for contig in contigs:
        header.add_meta('contig', items=[('ID', contig)])

    # TODO: add contig lengths?
    for info_items in VCF_HEADER_INFO_ITEMS:
        header.add_meta('INFO', items=[(k, str(v)) for k, v in info_items])

    return header


@contextmanager
def open_vcf_write(fp, contig: str) -> Generator[VariantFile, None, None]:
    with open_vcf(fp, mode='w', header=_get_vcf_header([contig])) as vcf:
        yield vcf
