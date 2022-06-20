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

from contextlib import closing, contextmanager
from functools import partial
import logging
import os
import shutil
from tempfile import NamedTemporaryFile
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from pysam import VariantFile, VariantHeader, VariantRecord
from ..enums import VariantType, VariantClassification
from ..utils import get_constant_category

VariantRecordTuple = Tuple[str, int, int, Optional[str], Optional[str], Optional[str], int, int]


VCF_MANIFEST_HEADER = {
    'vcf_alias',
    'vcf_id_tag',
    'vcf_path'
}

var_type_sub: int = VariantType.SUBSTITUTION.value
var_type_del: int = VariantType.DELETION.value
var_type_ins: int = VariantType.INSERTION.value
var_type_unk: int = VariantType.UNKNOWN.value

var_class_class: int = VariantClassification.CLASSIFIED.value
var_class_unclass: int = VariantClassification.UNCLASSIFIED.value
var_class_mono: int = VariantClassification.MONOMORPHIC.value


@contextmanager
def get_vcf(fp: str) -> Iterator[VariantFile]:
    with closing(VariantFile(fp)) as vcf:
        yield vcf


@contextmanager
def get_vcf_write(fp: str, header: VariantHeader) -> Iterator[VariantFile]:
    with closing(VariantFile(fp, mode='w', header=header)) as vcf:
        yield vcf


# TODO: add end_vcf for position 1 variants...?
def normalise_variant_record(
    get_var_id: Callable[[VariantRecord], Optional[str]],
    record: VariantRecord
) -> VariantRecordTuple:

    # Get record information
    chromosome: str = record.contig
    start: int = record.pos
    var_id: Optional[str] = get_var_id(record)
    ref: str = record.ref.upper()
    ref_len: int = len(ref)
    end: int = start + ref_len - 1

    if ref == '.':
        id_: str = f" {var_id}" if var_id is not None else ''
        raise ValueError(f"Empty REF for variant{id_} at {chromosome}:{start}!")

    # Monomorphic variant
    if record.alts is None:
        return chromosome, start, end, ref, None, var_id, var_type_unk, var_class_mono

    # Get alternative
    alt: str = record.alts[0].upper()
    alt_len: int = len(alt)

    if ref_len == 0 or alt_len == 0:
        raise ValueError("Invalid variant: REF and ALT must be set!")

    if ref_len == alt_len:
        return chromosome, start, end, ref, alt, var_id, var_type_sub, var_class_class
    else:
        # Check REF and ALT are preceded (or followed) by the same nucleotide
        pos_gt_one: bool = start > 1
        if ref_len > alt_len and alt_len == 1 and ref[0 if pos_gt_one else -1] == alt[0]:

            # Deletion
            if pos_gt_one:
                start += 1
                ref_trimmed = ref[1:]
            else:
                end -= 1
                ref_trimmed = ref[:-1]
            return chromosome, start, end, ref_trimmed, None, var_id, var_type_del, var_class_class

        elif ref_len < alt_len and ref_len == 1 and alt[0 if pos_gt_one else -1] == ref[0]:

            # Insertion
            if pos_gt_one:
                start += 1
                alt_trimmed = alt[1:]
            else:
                alt_trimmed = alt[:-1]
            end = start
            return chromosome, start, end, None, alt_trimmed, var_id, var_type_ins, var_class_class

        else:

            # Unclassified (possibly indel)
            return chromosome, start, end, ref, alt, var_id, var_type_sub, var_class_unclass


def _filter_var(chromosomes: Dict[str, Tuple[int, int]], record: VariantRecord) -> bool:
    chromosome: str = record.contig
    if chromosome not in chromosomes:
        return False

    pos_min, pos_max = chromosomes[chromosome]
    return (
        pos_min <= record.pos <= pos_max
        and (
            record.alts is None
            or len(record.alts) == 1
        )
    )


# TODO: if indexed, fetch by chromosome?
def load_vcf(fp: str, chromosomes: Dict[str, Tuple[int, int]], vcf_id_tag: Optional[str] = None) -> pd.DataFrame:
    if not vcf_id_tag:
        logging.info("No INFO tag specified as variant identifier for VCF file '%s', falling back to the ID field." % fp)

    def get_var_id_from_info(record: VariantRecord) -> Optional[str]:
        vcf_var_id: Any = record.info.get(vcf_id_tag, None)
        return str(vcf_var_id) if vcf_var_id is not None else None

    def get_var_id_from_id(record: VariantRecord) -> Optional[str]:
        return str(record.id) if record.id is not None else None

    get_var: Callable[[VariantRecord], VariantRecordTuple] = partial(
        normalise_variant_record,
        get_var_id_from_info if vcf_id_tag else get_var_id_from_id)

    filter_var: Callable[[VariantRecord], bool] = partial(
        _filter_var, chromosomes)

    with get_vcf(fp) as vcf:
        if vcf_id_tag and vcf_id_tag not in vcf.header.info:
            raise ValueError(f"INFO tag '{vcf_id_tag}' not found in VCF file '{fp}'!")

        df: pd.DataFrame = pd.DataFrame.from_records([
            get_var(record)
            for record in vcf.fetch()
            if filter_var(record)
        ], columns=[
            'Chromosome',
            'Start',
            'End',
            'ref',
            'alt',
            'vcf_var_id',
            'var_type',
            'var_class'
        ])

    # Compress information
    df.Chromosome = df.Chromosome.astype('category')
    df.Start = df.Start.astype(np.int32)
    df.End = df.End.astype(np.int32)
    df.ref = df.ref.astype('category')
    df.alt = df.alt.astype('category')
    df.vcf_var_id = df.vcf_var_id.astype('string')
    df.var_type = df.var_type.astype(np.int8)
    df.var_class = df.var_class.astype(np.int8)

    return df


def load_vcf_manifest(fp: str, chromosomes: Dict[str, Tuple[int, int]]) -> pd.DataFrame:

    def load_vcf_manifest_item(r: pd.Series) -> pd.DataFrame:
        vcf_id_tag: Optional[str] = r.vcf_id_tag if not pd.isnull(r.vcf_id_tag) else None
        df = load_vcf(r.vcf_path, chromosomes, vcf_id_tag=vcf_id_tag)
        df['vcf_alias'] = get_constant_category(
            r.vcf_alias, df.shape[0], categories=vcf_aliases)
        return df

    manifest_df: pd.DataFrame = pd.read_csv(fp)

    if set(manifest_df.columns) != VCF_MANIFEST_HEADER:
        raise ValueError("Invalid header!")

    if manifest_df.vcf_alias.duplicated(keep=False).any():
        raise ValueError("Duplicate VCF file aliases in manifest!")

    if manifest_df.vcf_path.duplicated(keep=False).any():
        raise ValueError("Duplicate VCF file paths in manifest!")

    logging.debug(f"{manifest_df.shape[0]} VCF files listed in manifest.")

    for fp in manifest_df.vcf_path.to_numpy():
        if not os.path.isfile(fp):
            raise FileNotFoundError(f"VCF file not found: '{fp}'!")

    vcf_aliases: List[str] = list(manifest_df.vcf_alias.to_numpy())
    return pd.concat(
        list(manifest_df.apply(load_vcf_manifest_item, axis=1)),
        ignore_index=True)


def _get_vcf_header_info_items(id_: str, type_: str, n: int) -> List[Tuple[str, Union[str, int]]]:
    return [
        ('ID', id_),
        ('Type', type_),
        ('Number', n)
    ]


VCF_HEADER_INFO_ITEMS: List[List[Tuple[str, Union[str, int]]]] = [
    _get_vcf_header_info_items(id_, type_, n)
    for id_, type_, n in [
        ('SGE_SRC', 'String', 1),
        ('SGE_REF', 'String', 1),
        ('SGE_OLIGO', 'String', 1),
        ('SGE_VCF_ALIAS', 'String', 1),
        ('SGE_VCF_VAR_ID', 'String', 1)
    ]
]


def _get_vcf_header(contigs: List[str]) -> VariantHeader:
    header: VariantHeader = VariantHeader()
    header.add_meta('source', value='valiant')

    for contig in contigs:
        header.add_meta('contig', items=[('ID', contig)])

    # TODO: add contig lengths?
    for info_items in VCF_HEADER_INFO_ITEMS:
        header.add_meta('INFO', items=info_items)

    return header


def write_vcf(fp: str, contigs: List[str], records: Iterable[Dict[str, Any]]) -> None:
    with NamedTemporaryFile('w', delete=True) as tmp:
        with get_vcf_write(tmp.name, _get_vcf_header(contigs)) as vcf:

            # Write records
            for record in records:
                vcf.write(vcf.new_record(**record))

        # Copy temporary file to final location
        shutil.copy(tmp.name, fp)
