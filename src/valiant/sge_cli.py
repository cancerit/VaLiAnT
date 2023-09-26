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

from itertools import chain
import logging
import os
import sys
from typing import Callable, Dict, Iterable, List, Optional, FrozenSet, Set, Tuple
import click
from pyranges import PyRanges

from .enums import TargetonMutator
from .errors import GenomicRangeOutOfBounds, InvalidBackground, InvalidVariantRef, SequenceNotFound
from .models.background_variants import GenomicPositionOffsets
from .models.base import GenomicRange, GenomicRangePair
from .models.codon_table import CodonTable
from .models.custom_variants import CustomVariant
from .models.dna_str import DnaStr
from .models.exon import CDSContextRepository, TranscriptInfo
from .models.annotation_repository import AnnotationRepository
from .models.exon_repository import ExonRepository
from .models.metadata_table import MetadataTable
from .models.new_pam import CdsPamBgAltSeqBuilder, PamBgAltSeqBuilder
from .models.oligo_generation_info import OligoGenerationInfo
from .models.oligo_segment import InvariantOligoSegment, OligoSegment, TargetonOligoSegment
from .models.oligo_template import OligoTemplate
from .models.options import Options
from .models.pam_protection import PamVariant
from .models.refseq_ranges import genomic_ranges_to_pyranges, ReferenceSequenceRangeCollection, ReferenceSequenceRanges, TargetReferenceRegion
from .models.refseq_repository import fetch_reference_sequences, ReferenceSequenceRepository
from .models.sequences import ReferenceSequence
from .models.sge_config import SGEConfig
from .models.snv_table import AuxiliaryTables
from .models.stats import Stats
from .models.targeton import ITargeton, PamProtCDSTargeton, PamProtTargeton
from .models.variant import BaseVariant, BaseVariantT
from .models.variant_repository_collection import VariantRepositoryCollection
from .common_cli import common_params, existing_file, finalise
from .cli_utils import load_codon_table
from .viz import get_text_diagram
from .writers import write_reference_sequences


def _load_gff_file(fp: Optional[str]) -> Optional[AnnotationRepository]:
    if fp:
        logging.debug("Loading GTF/GFF2 file...")
        try:
            return AnnotationRepository.from_gff(fp)
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to load GTF/GFF2 file!")
            sys.exit(1)
    else:
        logging.info("No GTF/GFF2 file provided.")
        return None


def _load_oligo_templates(fp: str) -> ReferenceSequenceRangeCollection:
    try:
        return ReferenceSequenceRangeCollection.load(fp)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical("Failed to load oligonucleotide templates!")
        sys.exit(1)


def resize_rsr(rsr: ReferenceSequenceRanges, gpo: GenomicPositionOffsets) -> ReferenceSequenceRanges:
    def resize_f(gr: GenomicRange) -> GenomicRange:
        return gr.resize(gpo.ref_to_alt_range(gr))

    return rsr.resize_regions(resize_f)


def get_oligo_template(
    rsr_ref: ReferenceSequenceRanges,
    ref: ReferenceSequenceRepository,
    cds: Optional[CDSContextRepository],
    bg_variants: Set[BaseVariant],
    pam_variants: Set[PamVariant],
    custom_variants: Set[CustomVariant],
    codon_table: CodonTable,
    config: SGEConfig
) -> OligoTemplate:

    bg_variant_ls = list(bg_variants)
    pam_variant_ls = list(pam_variants)

    # 1. Generate PAM protected reference sequence

    # Fetch reference sequence
    ref_seq = ref.get_genomic_range_sequence(rsr_ref.ref_range)
    if not ref_seq:
        raise ValueError("Sequence not found!")

    pam_ref_seq = PamBgAltSeqBuilder.from_ref(
        rsr_ref.ref_range, ref_seq, bg_variant_ls, pam_variant_ls)

    offsets: Optional[GenomicPositionOffsets] = None

    # TODO: to redesign the sequence splitting when ALT differs
    rsr = rsr_ref
    if bg_variants:
        offsets = GenomicPositionOffsets(
            pam_ref_seq.start,
            len(pam_ref_seq.ref_seq),
            pam_ref_seq.bg_variants
        )
        rsr_alt = resize_rsr(rsr_ref, offsets)
        rsr = rsr_alt

    # 2. Get constant regions

    def get_cds_targeton(
        region_pam_seq: PamBgAltSeqBuilder,
        exg_gr_pair: GenomicRangePair
    ) -> PamProtCDSTargeton:

        # Annotate region as CDS
        b = CdsPamBgAltSeqBuilder.from_noncds(
            region_pam_seq,
            get_cds_extension_sequence(exg_gr_pair[0]),
            get_cds_extension_sequence(exg_gr_pair[1]))

        # Validate background variants
        if not b.is_background_valid(
            codon_table,
            config.force_bg_fs,
            config.force_bg_ns
        ):
            raise InvalidBackground("Invalid background variants!")
        return PamProtCDSTargeton(b)

    def get_targeton(region_pam_seq: PamBgAltSeqBuilder) -> ITargeton:
        if cds:
            exg_gr_pair: Optional[GenomicRangePair] = cds.get_cds_extensions(
                region_pam_seq.pos_range)
            if exg_gr_pair is not None:
                return get_cds_targeton(region_pam_seq, exg_gr_pair)
        return PamProtTargeton(region_pam_seq)

    def get_constant_region_segment(gr: Optional[GenomicRange]) -> Optional[InvariantOligoSegment]:
        return InvariantOligoSegment(
            get_targeton(pam_ref_seq.get_pam_sub(gr))) if gr else None

    const_region_1: Optional[InvariantOligoSegment] = get_constant_region_segment(rsr.const_region_1)
    const_region_2: Optional[InvariantOligoSegment] = get_constant_region_segment(rsr.const_region_2)

    # 3. Get potential target regions

    def get_cds_extension_sequence(genomic_range: Optional[GenomicRange]) -> DnaStr:
        if genomic_range is None:
            return DnaStr.empty()

        is_local: bool = genomic_range.overlaps_range(
            pam_ref_seq.pos_range, unstranded=True)

        # TODO: should PAM variants be included?
        # Assumption: the local extension range is already in [local] cell-line coordinates
        sequence: Optional[str] = None
        try:
            sequence = (
                pam_ref_seq.get_pam_sub(genomic_range).pam_seq if is_local else
                ref.get_genomic_range_sequence(genomic_range)
            )
        except GenomicRangeOutOfBounds as ex:
            # TODO: handle CDS extension partially out of bounds
            logging.error(ex.args[0])

        if sequence is None:
            raise SequenceNotFound("CDS extension sequence not found!")

        return DnaStr(sequence)

    def get_oligo_segment(trr: TargetReferenceRegion) -> OligoSegment:
        region_pam_seq = pam_ref_seq.get_pam_sub(trr.genomic_range)
        return (
            TargetonOligoSegment(get_targeton(region_pam_seq), trr.mutators) if trr.mutators else
            get_constant_region_segment(trr.genomic_range)
        )

    def get_transcript_info(genomic_ranges: Iterable[GenomicRange]) -> Optional[TranscriptInfo]:
        if not cds:
            return None
        for genomic_range in genomic_ranges:
            transcript_info: Optional[TranscriptInfo] = cds.get_transcript_info(genomic_range)
            if transcript_info:
                return transcript_info
        return None

    # TODO: use ALT instead
    potential_target_segments: List[OligoSegment] = list(map(
        get_oligo_segment, rsr_ref.target_regions))

    # 4. Assemble regions in order
    segments: List[OligoSegment] = [const_region_1] if const_region_1 else []
    segments.extend(potential_target_segments)
    if const_region_2:
        segments.append(const_region_2)

    # Extract exon info (if any is available)
    transcript_info: Optional[TranscriptInfo] = get_transcript_info(
        segment.genomic_range for segment in segments)

    return OligoTemplate(
        rsr,
        offsets,
        transcript_info,
        pam_ref_seq,
        rsr.sgrna_ids,
        custom_variants,
        config.adaptor_5,
        config.adaptor_3,
        segments)


def get_oligo_templates(
    rsrs: ReferenceSequenceRangeCollection,
    ref: ReferenceSequenceRepository,
    variants: VariantRepositoryCollection,
    annotation: Optional[AnnotationRepository],
    codon_table: CodonTable,
    config: SGEConfig
) -> List[OligoTemplate]:

    def match_ref_regions_pam_variants() -> Dict[Tuple[str, int, int], Set[str]]:
        if not pam_variants_available:
            return {
                gr.as_unstranded(): set()
                for gr in rsrs.ref_ranges
            }

        # Map sgRNA variants to regions (ignoring their strandedness)
        range_variants: PyRanges = rsrs._ref_ranges.join(
            variants.pam._ranges, strandedness=False)

        # Count the number of variants per sgRNA per region for validation purposes
        # (chromosome, start, end) -> sgRNA ID -> variant count
        return {
            (chromosome, start + 1, end): set(rdf.sgrna_id.unique())
            for (chromosome, start, end), rdf in range_variants[[
                'sgrna_id',
                'variant_id'
            ]].as_df().groupby([
                'Chromosome',
                'Start',
                'End'
            ])
        }

    def get_missing_sgrna_err(genomic_range: GenomicRange, sgrna_ids: FrozenSet[str]) -> ValueError:
        return ValueError(
            "Expected sgRNA ID{} {} not found in region {}!".format(
                "'s" if len(sgrna_ids) > 1 else '',
                ', '.join(sgrna_ids),
                genomic_range.region)
        )

    def get_pam_variants(rsr: ReferenceSequenceRanges) -> Set[PamVariant]:
        if not rsr.sgrna_ids:
            return set()

        k: Tuple[str, int, int] = rsr.ref_range.as_unstranded()
        if k not in ref_ranges_sgrna_ids:
            raise get_missing_sgrna_err(rsr.ref_range, rsr.sgrna_ids)

        matching_sgrna_ids = ref_ranges_sgrna_ids[k]
        if rsr.sgrna_ids - matching_sgrna_ids:
            raise get_missing_sgrna_err(rsr.ref_range, rsr.sgrna_ids - matching_sgrna_ids)

        # Retrieve PAM protection variants
        return set(variants.pam.get_sgrna_variants_bulk(rsr.sgrna_ids))

    def get_bg_variants(rsr: ReferenceSequenceRanges) -> Set[BaseVariant]:
        return variants.background.get_variants(rsr.ref_range) if variants.background else set()

    def get_rsr_oligo_template(rsr: ReferenceSequenceRanges) -> OligoTemplate:

        def get_vars(label: str, cond: bool, f: Callable[[ReferenceSequenceRanges], Set[BaseVariantT]]) -> Set[BaseVariantT]:
            try:
                return f(rsr) if cond else set()
            except ValueError as ex:
                logging.critical(ex.args[0])
                logging.critical(
                    "Failed to retrieve %s variants "
                    "for region %s!" % (label, rsr.ref_range.region))
                sys.exit(1)

        bg_variants: Set[BaseVariant] = get_vars('background', bg_variants_available, get_bg_variants)
        pam_variants: Set[PamVariant] = get_vars('PAM protection', pam_variants_available, get_pam_variants)

        return get_oligo_template(
            rsr,
            ref,
            annotation.cds if annotation else None,
            bg_variants,
            pam_variants,
            variants.custom.get_variants(rsr.ref_range) if variants.custom else set(),
            codon_table,
            config)

    # Pass sgRNA variants to oligonucleotide template generation
    pam_variants_available: bool = variants.pam.count > 0
    bg_variants_available: bool = variants.background is not None
    ref_ranges_sgrna_ids = match_ref_regions_pam_variants()

    # Generate oligonucleotide templates
    ots: List[OligoTemplate] = list(map(get_rsr_oligo_template, rsrs._rsrs.values()))

    # Retrieve reference ranges of oligonucleotide templates with no transcript info (if any)
    ref_ranges_no_info: Set[GenomicRange] = set(
        ot.ref_range
        for ot in ots
        if not ot.transcript_info
    )

    # Return oligonucleotide templates if all have transcript information
    if not ref_ranges_no_info:
        return ots

    if annotation and annotation.utr:

        # Collect transcript information from UTR features
        utr_transcript_info: Dict[GenomicRange, TranscriptInfo] = annotation.utr.get_transcript_infos(
            genomic_ranges_to_pyranges(ref_ranges_no_info))

        if utr_transcript_info:

            utr_ranges: Set[GenomicRange] = set(utr_transcript_info.keys())
            ref_ranges_no_info -= utr_ranges

            # Log UTR matches
            for genomic_range in utr_ranges:
                logging.info(f"UTR found within region {genomic_range.region}.")

            # Set retrieved transcript information
            for ot in ots:
                if ot.ref_range in utr_ranges:
                    ot.transcript_info = utr_transcript_info[ot.ref_range]

    # Log missing transcript information if it could not be retrieved
    for genomic_range in ref_ranges_no_info:
        logging.info(f"No transcript information found for region {genomic_range.region}.")

    return ots


def get_oligo_template_qc_info(ot: OligoTemplate) -> List[str]:
    sequences: List[str] = [
        segment.sequence
        for segment in ot.ref_segments
    ]

    ranges: List[str] = [str(segment.start) for segment in ot.ref_segments]

    return [ot.name, ot.ref_range.region] + list(chain.from_iterable(zip(ranges, sequences)))


def generate_oligos(output: str, ref_repository: ReferenceSequenceRepository, aux: AuxiliaryTables, ot: OligoTemplate, species: str, assembly: str, options: Options) -> OligoGenerationInfo:

    # Generate metadata table
    mutation_df = ot.get_mutation_table(aux, options)
    metadata = MetadataTable.from_partial(
        species,
        assembly,
        mutation_df,
        options.oligo_min_length,
        options.oligo_max_length
    ) if mutation_df is not None else MetadataTable.empty()

    # Generate file name prefix
    base_fn = (
        '_'.join([ot.name] + sorted(ot.sgrna_ids)) if ot.sgrna_ids else
        ot.name
    )

    # Write to files
    metadata.write_sge_files(output, base_fn, ref_repository)

    # Log
    if metadata.short_oligo_n == 0:
        logging.warning(
            "Empty metadata table for targeton at %s (%s): no file generated!" % (
                ot.ref_range.region,
                ', '.join(ot.sgrna_ids) if ot.sgrna_ids else 'no PAM protection'
            ))
        logging.warning(
            "Empty unique oligonucleotides table for targeton at %s (%s): no file generated!" % (
                ot.ref_range.region,
                ', '.join(ot.sgrna_ids) if ot.sgrna_ids else 'no PAM protection'
            ))

    return metadata.get_info()


def run_sge(config: SGEConfig, sequences_only: bool) -> None:
    options = config.get_options()

    # Load codon table
    ct: CodonTable = load_codon_table(config.codon_table_fp)

    # Load CDS, stop codon, and UTR features from GTF/GFF2 file (if any)
    annotation: Optional[AnnotationRepository] = _load_gff_file(config.gff_fp)
    is_annotation_available: bool = annotation is not None

    # Load oligonucleotide templates
    rsrs: ReferenceSequenceRangeCollection = _load_oligo_templates(config.oligo_info_fp)

    # Collect all genomic ranges for which reference sequences have to be fetched
    ref_ranges = rsrs.ref_ranges

    variants = VariantRepositoryCollection.load(config, rsrs)

    if is_annotation_available and annotation.exons:

        # Retrieve CDS context (if any) for target regions
        # Only exonic sequences will have a CDS context
        try:
            annotation.cds = CDSContextRepository.from_exons(
                rsrs.target_ranges, annotation.exons)
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to match the CDS context!")
            sys.exit(1)
        except NotImplementedError as ex:
            logging.critical(ex.args[0])
            sys.exit(1)

        # Add CDS prefix and suffix ranges to ranges to fetch from reference
        ref_ranges |= annotation.cds.get_all_cds_extensions()

    # Initialise auxiliary tables
    all_mutators = rsrs.mutarors
    aux: AuxiliaryTables = AuxiliaryTables(
        ct,
        rsrs.strands,
        TargetonMutator.SNV in all_mutators,
        TargetonMutator.SNV_RE in all_mutators,
        TargetonMutator.AA in all_mutators)
    del all_mutators

    # Fetch reference sequences
    try:
        ref: ReferenceSequenceRepository = fetch_reference_sequences(config.ref_fasta_fp, ref_ranges)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical("Failed to retrieve reference sequences!")
        sys.exit(1)

    # Prepare oligonucleotides
    try:
        oligo_templates: List[OligoTemplate] = get_oligo_templates(
            rsrs,
            ref,
            variants,
            annotation,
            aux.codon_table,
            config)
    except (SequenceNotFound, InvalidBackground) as ex:
        logging.critical(ex.args[0])
        sys.exit(1)
    except InvalidVariantRef as ex:
        logging.critical(ex.args[0])
        if variants.background:
            logging.critical("Please check your background variants for overlaps!")
        sys.exit(1)

    def get_qc_row(ot: OligoTemplate) -> List[str]:
        return get_oligo_template_qc_info(ot)

    write_reference_sequences(
        map(get_qc_row, oligo_templates),
        os.path.join(config.output_dir, "ref_sequences.csv"))

    if sequences_only:
        sys.exit(0)

    stats = Stats()
    info: OligoGenerationInfo
    for ot in oligo_templates:
        try:

            # Generate all oligonucleotides and write to file
            info = generate_oligos(
                config.output_dir, ref, aux, ot, config.species, config.assembly, options)

        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to generate oligonucleotides!")
            sys.exit(1)

        stats.update(info)

    finalise(config, stats)


@click.command()
@common_params
@click.option('--gff', 'gff_fp', type=existing_file, help="Annotation GFF file path")
@click.option('--bg', 'bg_fp', type=existing_file, help="Background variant VCF file path")
@click.option('--pam', 'pam_fp', type=existing_file, help="PAM protection VCF file path")
@click.option('--vcf', 'vcf_fp', type=existing_file, help="Custom variant VCF manifest file path")
@click.option(
    '--sequences-only',
    is_flag=True,
    help="Fetch reference sequences and quit")
@click.option(
    '--revcomp-minus-strand',
    is_flag=True,
    help="Include reverse complement in oligonucleotide if reference is on minus strand")
@click.option(
    '--force-bg-ns',
    'force_bg_ns',
    is_flag=True,
    help="Allow non-synonymous background variants")
@click.option(
    '--force-bg-indels',
    'force_bg_fs',
    is_flag=True,
    help="Allow frame-shifting background variants")
def sge(

    # Input files
    oligo_info_fp: str,
    ref_fasta_fp: str,
    codon_table_fp: Optional[str],
    gff_fp: Optional[str],
    bg_fp: Optional[str],
    pam_fp: Optional[str],
    vcf_fp: Optional[str],

    # Output directory
    output_dir: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: Optional[str],
    adaptor_3: Optional[str],

    # Background variants
    force_bg_ns: bool,
    force_bg_fs: bool,

    # Actions
    sequences_only: bool,
    revcomp_minus_strand: bool,
    max_length: int,
    min_length: int

) -> None:
    """
    SGE oligonucleotide generation tool

    \b
    OLIGO_INFO is the BED-like input file path
    REF_FASTA is the reference genome FASTA file path
    OUTPUT is the output directory path
    SPECIES will be included in the metadata
    ASSEMBLY will be included in the metadata
    """

    run_sge(
        SGEConfig(
            species=species,
            assembly=assembly,
            adaptor_5=adaptor_5,
            adaptor_3=adaptor_3,
            min_length=min_length,
            max_length=max_length,
            codon_table_fp=codon_table_fp,
            oligo_info_fp=oligo_info_fp,
            ref_fasta_fp=ref_fasta_fp,
            output_dir=output_dir,
            revcomp_minus_strand=revcomp_minus_strand,
            force_bg_ns=force_bg_ns,
            force_bg_fs=force_bg_fs,
            gff_fp=gff_fp,
            bg_fp=bg_fp,
            pam_fp=pam_fp,
            vcf_fp=vcf_fp),
        sequences_only)
