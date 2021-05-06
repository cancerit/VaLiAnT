########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
#
# The Software remains the property of Genome Research Ltd ("GRL").
#
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
#
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
#
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
#
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
#############################

from itertools import chain
import logging
import os
import sys
from typing import Dict, Iterable, List, Optional, FrozenSet, Set, Tuple
import click
import pandas as pd
from pyranges import PyRanges
from .constants import DEFAULT_CODON_TABLE_FILE_NAME, DEFAULT_OLIGO_MAX_LENGTH, METADATA_FIELDS, METADATA_FIELDS_SET
from .enums import TargetonMutator
from .loaders.vcf import write_vcf
from .models.base import GenomicRange
from .models.codon_table import CodonTable
from .models.exon import AnnotationRepository, CDSContextRepository, GenomicRangePair, TranscriptInfo
from .models.metadata_table import MetadataTable
from .models.oligo_generation_info import OligoGenerationInfo
from .models.oligo_template import InvariantOligoSegment, OligoSegment, OligoTemplate, TargetonOligoSegment
from .models.options import Options
from .models.pam_protection import compute_pam_protected_sequence, PamProtectedReferenceSequence, PamProtectionVariantRepository, PamVariant
from .models.refseq_ranges import genomic_ranges_to_pyranges, ReferenceSequenceRangeCollection, ReferenceSequenceRanges, TargetReferenceRegion
from .models.refseq_repository import fetch_reference_sequences, ReferenceSequenceRepository
from .models.sequences import ReferenceSequence
from .models.snv_table import AuxiliaryTables
from .models.targeton import BaseTargeton, CDSTargeton, Targeton
from .models.variant import CustomVariant, VariantRepository
from .common_cli import common_params, existing_file
from .cli_utils import load_codon_table, validate_adaptor, set_logger
from .utils import get_constant_category, get_data_file_path, is_dna
from .writers import write_oligo_metadata, write_oligo_unique, write_reference_sequences


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


def _load_pam_protection_vcf(sgrna_ids: FrozenSet[str], pam: Optional[str]) -> PamProtectionVariantRepository:
    logging.debug("sgRNA ID's: %s" % ', '.join(sorted(sgrna_ids)) if sgrna_ids else "No sgRNA ID's.")

    vr = PamProtectionVariantRepository(sgrna_ids=set(sgrna_ids))

    if pam:
        if len(sgrna_ids) > 0:
            try:
                vr.load(pam)
            except ValueError as ex:
                logging.critical(ex.args[0])
                logging.critical("Failed to load the PAM protection variants!")
                sys.exit(1)
        else:
            logging.warning(
                "PAM protection variant file ignored: "
                "no sgRNA ID's associated to any targeton!")
    return vr


def get_oligo_template(
    rsr: ReferenceSequenceRanges,
    ref: ReferenceSequenceRepository,
    cds: Optional[CDSContextRepository],
    pam_variants: Set[PamVariant],
    custom_variants: Set[CustomVariant],
    adaptor_5: Optional[str] = None,
    adaptor_3: Optional[str] = None
) -> OligoTemplate:

    # 1. Generate PAM protected reference sequence

    def get_sequence(genomic_range: GenomicRange) -> ReferenceSequence:
        seq: Optional[str] = ref.get_genomic_range_subsequence(
            rsr.ref_range, genomic_range.start, genomic_range.end)
        if not seq:
            raise ValueError("Sequence not found!")
        return ReferenceSequence(seq, genomic_range)

    pam_ref_seq: PamProtectedReferenceSequence = compute_pam_protected_sequence(
        get_sequence(rsr.ref_range), pam_variants)

    # 2. Get constant regions

    def get_constant_segment(gr: Optional[GenomicRange]) -> Optional[InvariantOligoSegment]:
        return InvariantOligoSegment(pam_ref_seq.get_subsequence(gr)) if gr else None

    const_region_1: Optional[InvariantOligoSegment] = get_constant_segment(rsr.const_region_1)
    const_region_2: Optional[InvariantOligoSegment] = get_constant_segment(rsr.const_region_2)

    # 3. Get potential target regions

    def get_cds_extension_sequence(genomic_range: Optional[GenomicRange]) -> str:
        if genomic_range is None:
            return ''

        sequence: Optional[str] = ref.get_genomic_range_sequence(genomic_range)
        if sequence is None:
            raise RuntimeError("CDS extension sequence not found!")

        return sequence

    def get_cds_targeton(
        region_pam_seq: PamProtectedReferenceSequence,
        exg_gr_pair: GenomicRangePair
    ) -> CDSTargeton:
        return CDSTargeton.from_pam_seq(
            region_pam_seq,
            get_cds_extension_sequence(exg_gr_pair[0]),
            get_cds_extension_sequence(exg_gr_pair[1]))

    def get_targeton(region_pam_seq: PamProtectedReferenceSequence) -> BaseTargeton:
        if cds:
            exg_gr_pair: Optional[GenomicRangePair] = cds.get_cds_extensions(
                region_pam_seq.genomic_range)
            if exg_gr_pair is not None:
                return get_cds_targeton(region_pam_seq, exg_gr_pair)
        return Targeton.from_pam_seq(region_pam_seq)

    def get_oligo_segment(trr: TargetReferenceRegion) -> OligoSegment:
        region_pam_seq = pam_ref_seq.get_subsequence(trr.genomic_range)
        return (
            TargetonOligoSegment(get_targeton(region_pam_seq), trr.mutators) if trr.mutators else
            get_constant_segment(trr.genomic_range)
        )

    def get_transcript_info(genomic_ranges: Iterable[GenomicRange]) -> Optional[TranscriptInfo]:
        if not cds:
            return None
        for genomic_range in genomic_ranges:
            transcript_info: Optional[TranscriptInfo] = cds.get_transcript_info(genomic_range)
            if transcript_info:
                return transcript_info
        return None

    potential_target_segments: List[OligoSegment] = list(map(
        get_oligo_segment, rsr.target_regions))

    # 4. Assemble regions in order
    segments: List[OligoSegment] = [const_region_1] if const_region_1 else []
    segments.extend(potential_target_segments)
    if const_region_2:
        segments.append(const_region_2)

    # Extract exon info (if any is available)
    transcript_info: Optional[TranscriptInfo] = get_transcript_info(
        segment.genomic_range for segment in segments)

    return OligoTemplate(
        transcript_info,
        pam_ref_seq,
        rsr.sgrna_ids,
        custom_variants,
        adaptor_5,
        adaptor_3,
        segments)


def get_oligo_templates(
    rsrs: ReferenceSequenceRangeCollection,
    ref: ReferenceSequenceRepository,
    pam: PamProtectionVariantRepository,
    variant_repository: Optional[VariantRepository],
    annotation: Optional[AnnotationRepository],
    adaptor_5: Optional[str] = None,
    adaptor_3: Optional[str] = None
) -> List[OligoTemplate]:

    def match_ref_regions_pam_variants() -> Dict[Tuple[str, int, int], Set[str]]:
        if not pam_variants_available:
            return {
                gr.as_unstranded(): set()
                for gr in rsrs.ref_ranges
            }

        # Map sgRNA variants to regions (ignoring their strandedness)
        range_variants: PyRanges = rsrs._ref_ranges.join(pam._ranges, strandedness=False)

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
        return set(pam.get_sgrna_variants_bulk(rsr.sgrna_ids))

    def get_rsr_oligo_template(rsr: ReferenceSequenceRanges) -> OligoTemplate:
        try:
            pam_variants: Set[PamVariant] = get_pam_variants(rsr) if pam_variants_available else set()
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical(
                "Failed to retrieve PAM protection variants "
                "for region %s!" % rsr.ref_range.region)
            sys.exit(1)

        return get_oligo_template(
            rsr,
            ref,
            annotation.cds if annotation else None,
            pam_variants,
            variant_repository.get_variants(rsr.ref_range) if variant_repository else set(),
            adaptor_5=adaptor_5,
            adaptor_3=adaptor_3)

    # Pass sgRNA variants to oligonucleotide template generation
    pam_variants_available: bool = pam.count > 0
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


def get_oligo_template_qc_info(
    ref: ReferenceSequenceRepository,
    ot: OligoTemplate
) -> List[str]:

    sequences: List[str] = [
        segment.sequence
        for segment in ot.ref_segments
    ]

    ranges: List[str] = [str(segment.start) for segment in ot.ref_segments]

    if ref.get_genomic_range_sequence(ot.ref_range) != ''.join(sequences):
        logging.critical("Invalid reference sequence splitting!")
        sys.exit(1)

    return [ot.name, ot.ref_range.region] + list(chain.from_iterable(zip(ranges, sequences)))


def generate_oligos(output: str, ref_repository: ReferenceSequenceRepository, aux: AuxiliaryTables, ot: OligoTemplate, species: str, assembly: str, options: Options) -> OligoGenerationInfo:

    # Generate metadata table
    metadata = MetadataTable.from_partial(
        species,
        assembly,
        ot.get_mutation_table(aux, options),
        options.oligo_max_length)

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


@click.command()
@common_params
@click.option('--gff', type=existing_file)
@click.option('--pam', type=existing_file)
@click.option('--vcf', type=existing_file)
@click.option(
    '--sequences-only',
    is_flag=True,
    help="Fetch reference sequences and quit")
@click.option(
    '--revcomp-minus-strand',
    is_flag=True,
    help="Include reverse complement in oligonucleotide if reference is on minus strand")
def sge(

    # Input files
    oligo_info: str,
    ref_fasta: str,
    codon_table: Optional[str],
    gff: Optional[str],
    pam: Optional[str],
    vcf: Optional[str],

    # Output directory
    output: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: Optional[str],
    adaptor_3: Optional[str],

    # Actions
    sequences_only: bool,
    revcomp_minus_strand: bool,
    max_length: int,

    # Extra
    log: str

) -> None:
    """
    SGE oligonucleotide generation tool

    OLIGO_INFO is the BED-like input file path
    REF_FASTA is the reference genome FASTA file path
    OUTPUT is the output directory path
    """

    options: Options = Options(revcomp_minus_strand, max_length)

    # Set logging up
    set_logger(log)

    # Validate adaptor sequences
    validate_adaptor(adaptor_5)
    validate_adaptor(adaptor_3)

    # Load codon table
    ct: CodonTable = load_codon_table(codon_table)

    # Load CDS, stop codon, and UTR features from GTF/GFF2 file (if any)
    annotation: Optional[AnnotationRepository] = _load_gff_file(gff)
    exons: Optional[CDSContextRepository] = annotation.cds if annotation else None

    # Load oligonucleotide templates
    rsrs: ReferenceSequenceRangeCollection = _load_oligo_templates(oligo_info)

    # Load PAM protection variants
    pam_repository: PamProtectionVariantRepository = _load_pam_protection_vcf(rsrs.sgrna_ids, pam)

    # Collect all genomic ranges for which reference sequences have to be fetched
    ref_ranges = rsrs.ref_ranges

    # Load custom variants
    variant_repository: Optional[VariantRepository] = None
    if vcf:
        logging.debug("Loading custom variants...")
        try:
            variant_repository = VariantRepository.load(vcf, rsrs._ref_ranges)
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to load custom variants!")
            sys.exit(1)
        except FileNotFoundError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to load custom variants!")
            sys.exit(1)

    if exons:
        exons.register_target_ranges(rsrs.target_ranges)

        # Retrieve CDS context (if any) for target regions
        # Only exonic sequences will have a CDS context
        try:
            exons.compute_cds_contexts()
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to match the CDS context!")
            sys.exit(1)
        except NotImplementedError as ex:
            logging.critical(ex.args[0])
            sys.exit(1)

        # Add CDS prefix and suffix ranges to ranges to fetch from reference
        ref_ranges |= exons.get_all_cds_extensions()

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
        ref: ReferenceSequenceRepository = fetch_reference_sequences(ref_fasta, ref_ranges)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical("Failed to retrieve reference sequences!")
        sys.exit(1)

    # Prepare oligonucleotides
    oligo_templates: List[OligoTemplate] = get_oligo_templates(
        rsrs,
        ref,
        pam_repository,
        variant_repository,
        annotation,
        adaptor_5=adaptor_5,
        adaptor_3=adaptor_3)

    def get_qc_row(ot: OligoTemplate) -> List[str]:
        return get_oligo_template_qc_info(ref, ot)

    write_reference_sequences(
        map(get_qc_row, oligo_templates),
        os.path.join(output, "ref_sequences.csv"))

    if sequences_only:
        sys.exit(0)

    # Long oligonucleotides counter
    long_oligo_n: int = 0

    for ot in oligo_templates:
        try:
            # Generate all oligonucleotides and write to file
            info = generate_oligos(output, ref, aux, ot, species, assembly, options)
        except ValueError as ex:
            logging.critical(ex.args[0])
            logging.critical("Failed to generate oligonucleotides!")
            sys.exit(1)

        long_oligo_n += info.long_oligo_n

    # Log number of oligonucleotides discarded due to excessive length
    if long_oligo_n:
        logging.warning(
            "%d oligonucleotides longer than %d bases were discarded!" %
            (long_oligo_n, options.oligo_max_length))
