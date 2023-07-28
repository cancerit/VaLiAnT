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
from collections import namedtuple
from dataclasses import dataclass
from functools import partial
import logging
from typing import Dict, Iterable, List, Optional, Set, Tuple, ClassVar, Any, Callable, TypeVar
import pandas as pd
from pyranges import PyRanges
from pysam import VariantRecord

from valiant.constants import META_NEW, META_PAM_CODON_ALT, META_PAM_CODON_REF, META_REF, META_VCF_ALIAS, META_VCF_VAR_ID, METADATA_PAM_FIELDS
from .base import GenomicPosition, GenomicRange, PositionRange, StrandedPositionRange
from .refseq_repository import ReferenceSequenceRepository
from .sequences import ReferenceSequence
from ..enums import VariantType
from ..loaders.vcf import load_vcf, load_vcf_manifest, var_type_sub, var_type_del, var_type_ins, var_class_unclass, var_class_mono
from ..string_mutators import delete_nucleotides, insert_nucleotides, replace_nucleotides
from ..utils import get_id_column, is_dna, get_var_types, opt_str_length

# Metadata table fields used to generate the VCF output
VCF_RECORD_METADATA_FIELDS: List[str] = [
    'oligo_name',
    'ref',
    'new',
    'ref_seq',
    'pam_seq',
    'mut_position',
    'ref_start',
    'mutator',
    'ref_chr',
    'vcf_alias',
    'vcf_var_id',
    'var_type',
    *METADATA_PAM_FIELDS
]


T = TypeVar('T')


def _validate_seq(seq: str, label: str) -> None:
    if not seq:
        raise ValueError(f"Empty {label} sequence!")
    if not is_dna(seq):
        raise ValueError(f"Invalid {label} sequence '{seq}'!")


def _validate_ref(ref: str) -> None:
    _validate_seq(ref, 'reference')


def _validate_alt(alt: str) -> None:
    _validate_seq(alt, 'alternative')


def _validate_ref_in_target(seq: str, offset: int, ref: str) -> None:
    if seq[offset:offset + len(ref)] != ref:
        raise RuntimeError(f"Invalid variant: expected {ref}, found {seq[offset:offset + len(ref)]}!")


@dataclass(frozen=True)
class BaseVariant:
    __slots__ = {'genomic_position'}

    genomic_position: GenomicPosition

    type: ClassVar[VariantType]

    def get_ref(self) -> Optional[str]:
        return getattr(self, 'ref', None)

    def get_alt(self) -> Optional[str]:
        return getattr(self, 'alt', None)

    @property
    def ref_length(self) -> int:
        return opt_str_length(self.get_ref())

    @property
    def alt_length(self) -> int:
        return opt_str_length(self.get_alt())

    @property
    def alt_ref_delta(self) -> int:
        return self.alt_length - self.ref_length

    @property
    def is_frame_shifting(self) -> bool:
        return self.alt_ref_delta % 3 != 0

    @property
    def start(self) -> int:
        return self.genomic_position.position

    @property
    def ref_end(self) -> int:
        return self.start + min(0, self.ref_length - 1)

    @property
    def ref_start_pos(self) -> GenomicPosition:
        return self.genomic_position

    @property
    def ref_end_pos(self) -> GenomicPosition:
        return GenomicPosition(self.genomic_position.chromosome, self.ref_end)

    @property
    def ref_range(self) -> PositionRange:
        return PositionRange(self.start, self.ref_end)

    def in_range(self, genomic_range: GenomicRange) -> bool:
        return (
            self.ref_start_pos.in_range(genomic_range) or
            self.ref_end_pos.in_range(genomic_range)
        )

    def get_ref_offset(self, ref_seq: ReferenceSequence) -> int:
        if not ref_seq.genomic_range.contains_position(self.genomic_position):
            raise ValueError(
                f"Variant at {self.genomic_position} "
                f"not in genomic range {ref_seq.genomic_range.region}!")

        return self.genomic_position.position - ref_seq.genomic_range.start

    def _get_relative_position(self, seq_start: int) -> int:
        if seq_start < 1:
            raise ValueError("Invalid sequence start position!")
        return self.genomic_position.position - seq_start

    def get_corrected_ref_from(self, seq: str, offset: int) -> Optional[str]:
        return (
            seq[offset] if self.ref_length == 1 else
            seq[offset:offset + self.ref_length] if self.ref_length > 1 else
            None
        )

    def get_corrected_ref(self, seq: str, seq_start: int) -> Optional[str]:
        return self.get_corrected_ref_from(
            seq, self._get_relative_position(seq_start))

    def mutate_from(self, seq: str, offset: int, ref_check: bool = False) -> str:
        raise NotImplementedError()

    def mutate(self, seq: str, seq_start: int, ref_check: bool = False) -> str:
        return self.mutate_from(
            seq, self._get_relative_position(seq_start), ref_check=ref_check)


BaseVariantT = TypeVar('BaseVariantT', bound='BaseVariant')


def sort_variants(variants: Iterable[BaseVariantT]) -> List[BaseVariantT]:
    """Sort variants by genomic position"""

    return sorted(variants, key=lambda x: x.genomic_position.position)


def apply_variants(ref_seq: ReferenceSequence, variants: List[BaseVariant], ref_check: bool = False) -> str:
    alt_seq: str = ref_seq.sequence

    for variant in variants:

        # Validate variant genomic position relative to the sequence's
        offset: int = variant.get_ref_offset(ref_seq)

        # Update altered sequence
        alt_seq = variant.mutate_from(alt_seq, offset, ref_check=ref_check)

    return alt_seq


@dataclass(frozen=True)
class SubstitutionVariant(BaseVariant):
    __slots__ = {'genomic_position', 'ref', 'alt'}

    ref: str
    alt: str

    type: ClassVar[VariantType] = VariantType.SUBSTITUTION

    def __post_init__(self) -> None:
        _validate_ref(self.ref)
        _validate_alt(self.alt)

    @classmethod
    def from_variant_record(cls, r: VariantRecord) -> SubstitutionVariant:
        if not r.ref or not r.alts or not r.alts[0]:
            raise ValueError("Not a substitution!")
        return cls(GenomicPosition(r.contig, r.pos), r.ref, r.alts[0])

    def get_pyrange_record(self) -> Tuple[str, int, int]:
        position: int = self.genomic_position.position - 1
        return self.genomic_position.chromosome, position, position

    def mutate_from(self, seq: str, offset: int, ref_check: bool = False) -> str:
        if ref_check:
            _validate_ref_in_target(seq, offset, self.ref)
        return replace_nucleotides(seq, offset, self.ref, self.alt)


@dataclass(frozen=True)
class InsertionVariant(BaseVariant):
    __slots__ = {'genomic_position', 'alt'}

    alt: str

    type: ClassVar[VariantType] = VariantType.INSERTION

    def __post_init__(self) -> None:
        _validate_alt(self.alt)

    def mutate_from(self, seq: str, offset: int, ref_check: bool = False) -> str:
        return insert_nucleotides(seq, offset, self.alt)


@dataclass(frozen=True)
class DeletionVariant(BaseVariant):
    __slots__ = {'genomic_position', 'ref'}

    ref: str

    type: ClassVar[VariantType] = VariantType.DELETION

    def __post_init__(self) -> None:
        _validate_ref(self.ref)

    def mutate_from(self, seq: str, offset: int, ref_check: bool = False) -> str:
        if ref_check:
            _validate_ref_in_target(seq, offset, self.ref)
        return delete_nucleotides(seq, offset, self.ref)


@dataclass(frozen=True)
class CustomVariant:
    __slots__ = {'base_variant', 'vcf_alias', 'vcf_variant_id'}

    base_variant: BaseVariant
    vcf_alias: Optional[str]
    vcf_variant_id: Optional[str]

    def get_ref(self) -> Optional[str]:
        return self.base_variant.get_ref()

    @property
    def ref_length(self) -> int:
        return self.base_variant.ref_length

    def get_ref_range(self, strand: str) -> StrandedPositionRange:
        start: int = self.base_variant.genomic_position.position
        ref_length: int = len(getattr(self.base_variant, 'ref', ''))
        end: int = start + (ref_length if ref_length > 1 else 0)
        return StrandedPositionRange(start, end, strand)


VAR_TYPE_CONSTRUCTOR: Dict[int, Callable[[Any], BaseVariant]] = {
    var_type_sub: lambda t: SubstitutionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.ref, t.alt),
    var_type_del: lambda t: DeletionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.ref),
    var_type_ins: lambda t: InsertionVariant(
        GenomicPosition(t.Chromosome, t.Start_var + 1), t.alt)
}


def _map_variants(variants: pd.DataFrame, var_type: int, mask: bool = True) -> Dict[int, CustomVariant]:
    constructor = VAR_TYPE_CONSTRUCTOR[var_type]
    return {
        t.variant_id: CustomVariant(
            constructor(t),
            t.vcf_alias,
            t.vcf_var_id)
        for t in (
            variants[variants.var_type == var_type] if mask else
            variants
        ).itertuples(index=False)
    }


def _regions_to_chromosome_boundaries(regions: PyRanges) -> Dict[str, Tuple[int, int]]:
    return {
        chromosome: (df.Start.min() + 1, df.End.max())
        for chromosome, df in regions.dfs.items()
    }


@dataclass
class VariantRepository:
    _variants: Dict[int, CustomVariant]
    _region_variants: Dict[Tuple[str, int, int], Set[int]]

    @classmethod
    def load_vcf(cls, vcf_fp: str, regions: PyRanges) -> VariantRepository:
        # TODO: verify variant filtering criteria being applied at this stage (if any)
        chromosome_boundaries = _regions_to_chromosome_boundaries(regions)
        df = load_vcf(vcf_fp, chromosome_boundaries)
        # TODO: discuss variant identifiers
        df['vcf_alias'] = None
        df['vcf_var_id'] = None
        return cls.from_df(regions, df)

    @classmethod
    def load(cls, manifest_fp: str, regions: PyRanges) -> VariantRepository:

        # Load all permitted variants from multiple VCF files
        chromosome_boundaries = _regions_to_chromosome_boundaries(regions)
        df = load_vcf_manifest(manifest_fp, chromosome_boundaries)
        return cls.from_df(regions, df)

    @classmethod
    def from_df(cls, regions: PyRanges, custom_variants: pd.DataFrame) -> VariantRepository:
        custom_variants_n: int = custom_variants.shape[0]

        if custom_variants_n == 0:
            return cls({}, {})

        logging.debug("Collected %d custom variants." % custom_variants_n)

        # Make start positions zero-based
        custom_variants.Start -= 1

        # Assign identifier to all variants
        custom_variants['variant_id'] = get_id_column(custom_variants_n)
        custom_variant_ranges: PyRanges = PyRanges(df=custom_variants)

        # Match regions with variants
        ref_ranges_variants: pd.DataFrame = regions.join(
            custom_variant_ranges,
            suffix='_var'
        )[[
            'Start_var',
            'End_var',
            'variant_id',
            'ref',
            'alt',
            'vcf_alias',
            'vcf_var_id',
            'var_type',
            'var_class'
        ]].as_df()
        del custom_variant_ranges

        variants: pd.DataFrame = ref_ranges_variants.drop([
            'Start',
            'End'
        ], axis=1).drop_duplicates([
            'variant_id'
        ])

        # Log and discard monomorphic variants
        mono_mask: pd.Series = variants.var_class == var_class_mono
        if mono_mask.any():
            for chromosome, start in variants.loc[mono_mask, [
                'Chromosome',
                'Start_var'
            ]].itertuples(index=False, name=None):
                logging.info(f"Monomorphic variant at {chromosome}:{start + 1} (SKIPPED).")
            variants = variants[~mono_mask]
            mono_mask = ref_ranges_variants.var_class != var_class_mono
            ref_ranges_variants = ref_ranges_variants[mono_mask]
        del mono_mask

        # Log unclassified variants
        unclass_mask: pd.Series = variants.var_class == var_class_unclass
        for r in variants[unclass_mask].itertuples(index=False):
            logging.info(f"Unclassified variant at {r.Chromosome}:{r.Start_var + 1}: {r.ref}>{r.alt}.")
        del unclass_mask

        # List all variant types represented in the collection
        var_types: List[int] = get_var_types(variants.var_type)
        var_types_n: int = len(var_types)

        if not var_types_n:
            raise RuntimeError("No variant types in variant table!")

        # Map variant indices to variant objects
        var_types_gt_one: bool = var_types_n > 1
        matching_variants: Dict[int, CustomVariant] = _map_variants(
            variants, var_types[0], mask=var_types_gt_one)
        if var_types_gt_one:
            for var_type in var_types[1:]:
                matching_variants.update(_map_variants(variants, var_type))
        del variants, var_types_gt_one

        # Map regions to variant indices
        ref_ranges_variant_ids: Dict[Tuple[str, int, int], Set[int]] = {
            (chromosome, start + 1, end): set(g.variant_id[(
                (g.Start_var >= start)
                & (g.End_var <= end)
            )].unique())
            for (chromosome, start, end), g in ref_ranges_variants.groupby([
                'Chromosome',
                'Start',
                'End'
            ])
        }

        return cls(matching_variants, ref_ranges_variant_ids)

    def get_variants(self, genomic_range: GenomicRange) -> Set[CustomVariant]:
        r: Tuple[str, int, int] = genomic_range.as_unstranded()

        if r not in self._region_variants or not self._region_variants[r]:
            return set()

        return set(self._variants[var_id] for var_id in self._region_variants[r])


def get_variant_from_tuple(chromosome: str, position: int, ref: str, alt: str) -> BaseVariant:
    genomic_position: GenomicPosition = GenomicPosition(chromosome, position)
    ref_len: int = len(ref)
    alt_len: int = len(alt)

    if ref_len == 0 or alt_len == 0:
        raise ValueError("Invalid variant: REF and ALT must be set!")

    if ref_len == alt_len:
        return SubstitutionVariant(genomic_position, ref, alt)
    else:
        # Check REF and ALT are preceded (or followed) by the same nucleotide
        pos_gt_one: bool = position > 1
        if ref_len > alt_len and alt_len == 1 and ref[0 if pos_gt_one else -1] == alt[0]:

            # Deletion
            ref_trimmed: str = ref[1:] if pos_gt_one else ref[:-1]
            if pos_gt_one:
                genomic_position += 1
            return DeletionVariant(genomic_position, ref_trimmed)

        elif ref_len < alt_len and ref_len == 1 and alt[0 if pos_gt_one else -1] == ref[0]:

            # Insertion
            alt_trimmed: str = alt[1:] if pos_gt_one else alt[:-1]
            if pos_gt_one:
                genomic_position += 1
            return InsertionVariant(genomic_position, alt_trimmed)

        else:

            # Unclassified (possibly indel)
            logging.info(f"Unclassified variant at {chromosome}:{position}: {ref}>{alt}.")
            return SubstitutionVariant(genomic_position, ref, alt)


def get_variant(r: VariantRecord) -> BaseVariant:
    return get_variant_from_tuple(r.contig, r.pos, r.ref, r.alts[0])


def _get_shared_nucleotide(
    ref_repository: ReferenceSequenceRepository,
    ref_seq: str,
    pam_seq: str,
    chromosome: str,
    seq_start: int,
    prev_pos: int
) -> Tuple[str, str]:
    offset: int = prev_pos - seq_start
    if offset >= 0:
        return pam_seq[offset], ref_seq[offset]
    else:
        shared_nt: str = ref_repository.get_nucleotide_unsafe(
            chromosome, prev_pos)
        return shared_nt, shared_nt


def _get_shared_nucleotide_and_slice(
    ref_repository: ReferenceSequenceRepository,
    ref_seq: str,
    pam_seq: str,
    chromosome: str,
    seq_start: int,
    prev_pos: int,
    mut_len: int
) -> Tuple[str, str, str]:
    shared_nt: str
    pam_slice: str
    ref_slice: str
    offset: int = prev_pos - seq_start

    if offset >= 0:

        # Get nucleotide shared between REF and ALT
        shared_nt = pam_seq[offset]

        # Retrieve reference sequence before and after PAM protection
        sl = slice(offset, offset + mut_len + 1)
        pam_slice = pam_seq[sl]
        ref_slice = ref_seq[sl]

    else:

        # Get nucleotide shared between REF and ALT
        shared_nt = ref_repository.get_nucleotide_unsafe(chromosome, prev_pos)

        # Retrieve reference sequence before and after PAM protection
        sl = slice(0, mut_len + 1)
        pam_slice = shared_nt + pam_seq[sl]
        ref_slice = shared_nt + ref_seq[sl]

    return shared_nt, ref_slice, pam_slice


def _get_offset_vcf_record(
    ref_repository: ReferenceSequenceRepository,
    chromosome: str,
    pos: int,
    seq_start: int,
    ref_seq: str,
    pam_seq: str,
    mut_len: int,
) -> Tuple[int, int, str, str, str]:
    start: int
    end: int
    ref_slice: str
    pam_slice: str
    shared_nt: str

    if pos == 1:
        start = 1
        end = mut_len + 2

        # Retrieve nucleotide shared between REF and ALT
        shared_nt = pam_seq[mut_len]

        # Retrieve reference sequence before and after PAM protection
        sl = slice(0, mut_len + 1)
        ref_slice = ref_seq[sl]
        pam_slice = pam_seq[sl]

    else:
        start = pos - 1
        end = start + mut_len + 1

        # Retrieve nucleotide shared between REF and ALT and reference sequence before PAM protection
        shared_nt, ref_slice, pam_slice = _get_shared_nucleotide_and_slice(
            ref_repository, ref_seq, pam_seq, chromosome, seq_start, start, mut_len)

    return start, end, shared_nt, ref_slice, pam_slice


def _get_insertion_record(
    ref_repository: ReferenceSequenceRepository,
    chromosome: str,
    pos: int,
    alt: str,
    seq_start: int,
    ref_seq: str,
    pam_seq: str
) -> Tuple[int, int, str, str, Optional[str]]:
    if pos == 1:
        pam_nt = pam_seq[0]
        ref_nt = ref_seq[0]
        alt_ = alt + pam_nt
        return 1, 2, pam_nt, alt_, (ref_nt if pam_nt != ref_nt else None)
    else:
        prev_pos: int = pos - 1
        pam_nt, ref_nt = _get_shared_nucleotide(
            ref_repository, ref_seq, pam_seq, chromosome, seq_start, prev_pos)
        alt_ = pam_nt + alt
        return prev_pos, pos, pam_nt, alt_, (ref_nt if pam_nt != ref_nt else None)


def _get_deletion_record(
    ref_repository: ReferenceSequenceRepository,
    chromosome: str,
    pos: int,
    ref: str,
    seq_start: int,
    ref_seq: str,
    pam_seq: str
) -> Tuple[int, int, str, str, Optional[str]]:
    mut_len: int = len(ref)
    prev_pos, end, shared_nt, ref_slice, pam_slice = _get_offset_vcf_record(
        ref_repository, chromosome, pos, seq_start, ref_seq, pam_seq, mut_len)
    return prev_pos, end, pam_slice, shared_nt, (ref_slice if pam_slice != ref_slice else None)


def _get_slices(ref_seq: str, pam_seq: str, offset: int, mut_len: int) -> Tuple[str, str]:
    if mut_len > 1:
        sl = slice(offset, offset + mut_len)
        ref_slice = ref_seq[sl]
        pam_slice = pam_seq[sl]
    else:
        ref_slice = ref_seq[offset]
        pam_slice = pam_seq[offset]
    return ref_slice, pam_slice


def _get_substitution_record(
    pos: int,
    ref: str,
    alt: str,
    seq_start: int,
    ref_seq: str,
    pam_seq: str
) -> Tuple[int, int, str, str, Optional[str]]:
    mut_len: int = len(ref)
    end: int = pos + mut_len

    # Retrieve reference sequence before and after PAM protection
    offset: int = pos - seq_start
    ref_slice, pam_slice = _get_slices(ref_seq, pam_seq, offset, mut_len)
    return pos, end, pam_slice, alt, (ref_slice if pam_slice != ref_slice else None)


def get_nullable_field(s: namedtuple, field: str, default: Optional[T] = None) -> Optional[T]:
    x: T = s.__getattribute__(field)
    return x if not pd.isna(x) else default


@dataclass
class VCFRecordParams:
    __slots__ = ['chr', 'pos', 'ref', 'alt', 'var_type', 'mutator', 'oligo_name', 'vcf_alias', 'vcf_var_id', 'ref_start', 'ref_seq', 'pam_seq']

    chr: str
    pos: int
    ref: Optional[str]
    alt: Optional[str]
    var_type: int
    mutator: str
    oligo_name: str
    vcf_alias: Optional[str]
    vcf_var_id: Optional[str]
    ref_start: int
    ref_seq: str
    pam_seq: str

    def __post_init__(self) -> None:
        assert not (self.ref is None and self.alt is None)
        assert isinstance(self.pos, int)
        assert isinstance(self.var_type, int)
        assert isinstance(self.ref_start, int)

    @classmethod
    def from_meta(cls, pam_mode: bool, meta: namedtuple) -> VCFRecordParams:
        f = (
            cls.from_meta_pam_ext if (
                pam_mode
                and meta.pam_codon_mask
                and meta.pam_mut_start is not None
            ) else
            cls.from_meta_no_pam
        )
        return f(meta)

    @classmethod
    def from_meta_no_pam(cls, meta: namedtuple) -> VCFRecordParams:
        return cls(
            chr=meta.ref_chr,
            pos=int(meta.mut_position),
            ref=get_nullable_field(meta, META_REF) or None,
            alt=get_nullable_field(meta, META_NEW) or None,
            var_type=int(meta.var_type),
            mutator=meta.mutator,
            oligo_name=meta.oligo_name,
            vcf_alias=get_nullable_field(meta, META_VCF_ALIAS),
            vcf_var_id=get_nullable_field(meta, META_VCF_VAR_ID),
            ref_start=int(meta.ref_start),
            ref_seq=meta.ref_seq,
            pam_seq=meta.pam_seq)

    @classmethod
    def from_meta_pam_ext(cls, meta: namedtuple) -> VCFRecordParams:
        ref = get_nullable_field(meta, META_PAM_CODON_REF) or None
        alt = get_nullable_field(meta, META_PAM_CODON_ALT) or None
        var_type = (
            var_type_del if alt is None else
            var_type_ins if ref is None else
            var_type_sub
        )
        return cls(
            chr=meta.ref_chr,
            pos=int(meta.pam_mut_start),
            ref=ref,
            alt=alt,
            var_type=int(var_type),
            mutator=meta.mutator,
            oligo_name=meta.oligo_name,
            vcf_alias=get_nullable_field(meta, META_VCF_ALIAS),
            vcf_var_id=get_nullable_field(meta, META_VCF_VAR_ID),
            ref_start=int(meta.ref_start),
            ref_seq=meta.ref_seq,
            pam_seq=meta.pam_seq)

    def get_vcf_record(
        self,
        ref_repository: ReferenceSequenceRepository,
        pam_mode: bool
    ) -> Dict[str, Any]:
        """
        Convert into a VCF record

        If pam_ref is true, use the PAM-protected sequence as reference,
        with extending to cover the whole start and end codons if those are
        PAM-protected; the original reference otherwise.
        """

        pos: int
        ref: str
        alt: str
        end: int
        pre_pam: Optional[str]  # Reference before PAM protection (only set if it differs from ref)

        # TODO: verify stop position for variants at position 1... might need correcting here as well
        if self.var_type == var_type_del:
            pos, end, ref, alt, pre_pam = _get_deletion_record(
                ref_repository, self.chr, self.pos, self.ref, self.ref_start, self.ref_seq, self.pam_seq)
        elif self.var_type == var_type_ins:
            pos, end, ref, alt, pre_pam = _get_insertion_record(
                ref_repository, self.chr, self.pos, self.alt, self.ref_start, self.ref_seq, self.pam_seq)
        elif self.var_type == var_type_sub:
            pos, end, ref, alt, pre_pam = _get_substitution_record(
                self.pos, self.ref, self.alt, self.ref_start, self.ref_seq, self.pam_seq)
        else:
            raise RuntimeError("Invalid variant type!")

        # Set INFO tags
        info: Dict[str, Any] = {
            'SGE_SRC': self.mutator,
            'SGE_OLIGO': self.oligo_name
        }
        if pam_mode and pre_pam:
            info['SGE_REF'] = pre_pam

        if self.vcf_var_id:
            info['SGE_VCF_ALIAS'] = self.vcf_alias
            info['SGE_VCF_VAR_ID'] = self.vcf_var_id

        # Set VCF record fields
        return {
            'alleles': (ref if pam_mode else (pre_pam or ref), alt),
            'contig': self.chr,
            'start': pos - 1,  # zero-based representation
            'stop': end - 1,  # zero-based representation
            'info': info
        }


def _get_vcf_record_params_constructor(pam_ref: bool) -> Callable:
    return (
        partial(VCFRecordParams.from_meta, pam_ref) if pam_ref else
        VCFRecordParams.from_meta_no_pam
    )


def get_records(
    ref_repository: ReferenceSequenceRepository,
    pam_ref: bool,
    meta: pd.DataFrame
) -> List[Dict[str, Any]]:

    to_record_params_f = _get_vcf_record_params_constructor(pam_ref)

    def get_record(x: namedtuple) -> Dict[str, Any]:
        return to_record_params_f(x).get_vcf_record(
            ref_repository, pam_ref)

    return list(map(get_record, meta[VCF_RECORD_METADATA_FIELDS].itertuples(index=False)))
