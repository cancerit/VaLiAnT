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

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator

from pysam import VariantFile, VariantHeader

from .loaders.vcf import open_vcf


@dataclass(slots=True)
class Allele:
    s: str
    nt: str | None = None

    def get_vcf_allele(self, pos: int) -> str:
        if not self.nt:
            return self.s
        return (
            # TODO: or zero?
            (self.nt + self.s) if pos > 1 else
            (self.s + self.nt)
        )

    def clone(self) -> Allele:
        return Allele(self.s, nt=self.nt)


@dataclass(slots=True)
class VcfVariant:
    start: int
    end: int
    ref: Allele
    alt: Allele
    sge_ref: Allele | None = None

    @classmethod
    def from_partial(
        cls,
        start: int,
        ref: Allele,
        alt: Allele,
        end: int | None = None,
        sge_ref: Allele | None = None
    ) -> VcfVariant:
        assert start >= 0
        assert end is None or end >= 0

        if end is None:
            end = start + max(0, len(ref.s) - 1)

        if not ref.s or not alt.s:
            # TODO: verify it may not be zero instead at this stage
            if start > 1:
                start -= 1
            else:
                # TODO: verify end correction
                end += 1

        return cls(start, end, ref, alt, sge_ref=sge_ref)

    def __post_init__(self) -> None:
        assert self.start >= 0
        assert self.end >= self.start
        assert self.ref
        assert self.alt


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
    with _open_vcf_write(fp, _get_vcf_header([contig])) as vcf:
        yield vcf


def get_vcf_record_info(
    variant: VcfVariant,
    oligo_name: str,
    mutator: str,
    vcf_alias: str | None,
    vcf_var_id: str | None
) -> dict[str, str]:
    vcf_info = {
        'SGE_SRC': mutator,
        'SGE_OLIGO': oligo_name
    }

    if variant.sge_ref:
        sge_ref = variant.sge_ref.get_vcf_allele(variant.start)
        ref = variant.ref.get_vcf_allele(variant.start)
        assert len(sge_ref) == len(ref)
        if sge_ref != ref:
            vcf_info['SGE_REF'] = sge_ref

    if vcf_alias and vcf_var_id:
        # TODO: this is for compatibility purposes (BUG in v. 3.*?)
        vcf_info['SGE_VCF_ALIAS'] = vcf_alias

    if vcf_var_id:
        vcf_info['SGE_VCF_VAR_ID'] = vcf_var_id

    return vcf_info


def write_vcf_record(
    vcf: VariantFile,
    contig: str,
    variant: VcfVariant,
    mutator: str,
    oligo_name: str,
    vcf_alias: str | None,
    vcf_var_id: str | None
) -> None:

    # Populate the INFO field
    vcf_info = get_vcf_record_info(
        variant,
        oligo_name,
        mutator,
        vcf_alias,
        vcf_var_id)

    ref = variant.ref.get_vcf_allele(variant.start)
    alt = variant.alt.get_vcf_allele(variant.start)

    # Write record to file
    vcf.write(vcf.new_record(
        alleles=(ref, alt),
        contig=contig,
        start=variant.start,
        # stop=variant.end,
        stop=variant.start + max(0, len(ref) - 1),
        info=vcf_info))
