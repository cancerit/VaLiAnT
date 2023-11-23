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
import os
from sqlite3 import Connection

import click
from pysam.libcfaidx import FastaFile

from .annot_variant import AnnotVariant
from .annotation import Annotation
from .codon_table import CodonTable
from .codon_table_loader import load_codon_table_rows
from .common_cli import common_params, existing_file
from .constants import OUTPUT_CONFIG_FILE_NAME
from .contig_filter import ContigFilter
from .custom_variant import CustomVariant
from .db import get_db_conn, cursor, dump_all
from .loaders.experiment import ExperimentConfig, TargetonConfig
from .loaders.fasta import open_fasta
from .loaders.gtf import GtfLoader
from .loaders.vcf_manifest import VcfManifest
from .mutator import MutatorCollection
from .pam_variant import PamVariant
from .pattern_variant import PatternVariant
from .queries import insert_annot_pattern_variants, insert_custom_variant_collection, insert_exons, insert_background_variants, \
    insert_pam_protection_edits, insert_gene_offsets, insert_pattern_variants, select_exons_in_range
from .seq import Seq
from .sge_config import SGEConfig
from .transcript_seq import TranscriptSeq
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .utils import get_default_codon_table_path, get_ddl_path


def fetch_sequence(fa: FastaFile, contig: str, r: UIntRange) -> Seq:
    # TODO: verify coordinate convention and need to extend for the purposes of VCF generation
    logging.debug("Fetching reference sequence at %s:%d-%d." %
                  (contig, r.start, r.end))
    return Seq(r.start, DnaStr(
        fa.fetch(reference=contig, start=r.start - 1, end=r.end).upper()))


def load_sequences(fa: FastaFile, contig: str, ranges: list[UIntRange]) -> list[Seq]:
    logging.debug("Fetching %d reference sequences..." % len(ranges))
    return [
        fetch_sequence(fa, contig, r)
        for r in ranges
    ]


def load_custom_variants(conn: Connection, vcf_manifest_fp: str, contig: str, targeton_ranges: list[UIntRange]) -> None:

    # Load manifest file
    manifest = VcfManifest.load(vcf_manifest_fp)
    manifest.test_vcf_extance()

    # Set positional filter
    contig_ft = ContigFilter.from_ranges(contig, targeton_ranges)

    total_custom_variants: int = 0
    for vcf_proxy in manifest.vcfs:

        # Load variants from VCF file
        custom_vars = CustomVariant.load_vcf(
            vcf_proxy.vcf_path,
            contig_ft,
            vcf_id_tag=vcf_proxy.vcf_id_tag)

        # Push the variants to the database
        insert_custom_variant_collection(
            conn, vcf_proxy.vcf_alias, custom_vars)

        total_custom_variants += len(custom_vars)

    logging.debug("Collected %d custom variants." % total_custom_variants)


def load_background_variants(conn: Connection, vcf_fp: str, contig: str, ranges: list[UIntRange]) -> None:

    # Set positional filter
    contig_ft = ContigFilter.from_ranges(contig, ranges)

    # Load variants from VCF file
    bg_vars = CustomVariant.load_vcf(vcf_fp, contig_ft)
    logging.debug("Collected %d background variants." % len(bg_vars))

    # Push the variants to the database
    insert_background_variants(conn, bg_vars)


def load_pam_variants(conn: Connection, vcf_fp: str, contig: str, ranges: list[UIntRange]) -> None:

    # Set positional filter
    contig_ft = ContigFilter.from_ranges(contig, ranges)

    # Load variants from VCF file
    pam_vars = PamVariant.load_vcf(vcf_fp, contig_ft)
    logging.debug("Collected %d PAM variants." % len(pam_vars))

    # Push the variants to the database
    insert_pam_protection_edits(conn, pam_vars)


def init_database(conn: Connection) -> None:
    with open(get_ddl_path()) as fh:
        with cursor(conn) as cur:
            cur.executescript(fh.read())


def get_targeton_region_exon_id(conn: Connection, r: UIntRange) -> int | None:
    exon_ids = select_exons_in_range(conn, r.start, r.end)
    en: int = len(exon_ids)
    match en:
        case 0:
            return None
        case 1:
            return exon_ids[0]
        case _:
            raise ValueError("Invalid targeton region: overlaps multiple exons!")


def get_pattern_variants_from_region(
    conn: Connection,
    codon_table: CodonTable,
    transcript: TranscriptSeq | None,
    targeton: Seq,
    region: UIntRange,
    mc: MutatorCollection
) -> tuple[list[PatternVariant], list[AnnotVariant]]:

    # Get overlapping exon ID's
    exon_id = get_targeton_region_exon_id(conn, region)
    is_cds = exon_id is not None

    if is_cds:
        assert transcript is not None
        r_seq = transcript.get_cds_seq(exon_id, region)
    else:
        r_seq = targeton.subseq(region, rel=False)

    vars, annot_vars = mc.get_variants(codon_table, r_seq)
    insert_pattern_variants(conn, vars)
    insert_annot_pattern_variants(conn, annot_vars)

    return vars, annot_vars


def process_targeton_config(
    conn: Connection,
    codon_table: CodonTable,
    transcript: TranscriptSeq | None,
    targetons: dict[int, Seq],
    t: TargetonConfig
) -> tuple[list[PatternVariant], list[AnnotVariant]]:
    pattern_variants: list[PatternVariant] = []
    annot_variants: list[AnnotVariant] = []

    # TODO: apply pattern variants to R1 and R3 as well
    for i, r in enumerate([None, t.region_2, None]):
        if r is not None:
            mutators = t.mutators[i]
            if mutators:
                targeton = targetons[t.ref.start]
                mc = MutatorCollection.from_configs(mutators)
                vars, annot_vars = get_pattern_variants_from_region(
                    conn, codon_table, transcript, targeton, r, mc)
                pattern_variants.extend(vars)
                annot_variants.extend(annot_vars)

    return pattern_variants, annot_variants


def get_background_context_range(targeton_ranges, annot: Annotation | None) -> UIntRange:
    """Get the minimum range in which to apply background variants"""

    # Identify minimal context required
    # TODO: extend to support multiple contigs
    # Assumption: single ContigFilter
    t_min = min(t.start for t in targeton_ranges)
    t_max = max(t.end for t in targeton_ranges)
    if annot:
        bg_ctx = UIntRange(min(t_min, annot.cds_start), max(t_max, annot.cds_end))
    else:
        bg_ctx = UIntRange(t_min, t_max)
    return bg_ctx


def run_sge(config: SGEConfig, sequences_only: bool) -> None:

    # Load options
    options = config.get_options()

    # Load codon table
    if not config.codon_table_fp:
        logging.info(
            "Codon table not specified, the default one will be used.")
    codon_table = CodonTable.from_list(load_codon_table_rows(
        config.codon_table_fp or get_default_codon_table_path()))

    # Load experiment configuration
    exp = ExperimentConfig.load(config.oligo_info_fp)
    targeton_ranges = exp.ref_ranges

    # Load annotation
    annot: Annotation | None = None
    if config.gff_fp:
        annot = GtfLoader(exp.contig, exp.strand).load_gtf(config.gff_fp)

    all_ranges = targeton_ranges + [
        UIntRange(x.start, x.end)
        for x in annot.cds.ranges
    ] if annot else []

    transcript: TranscriptSeq | None = None

    # Get FASTA file
    with open_fasta(config.ref_fasta_fp) as fa:
        # TODO: load all sequences at once
        targeton_ref_seqs = load_sequences(fa, exp.contig, targeton_ranges)
        if annot:
            cds_ref_seqs = load_sequences(fa, exp.contig, annot.cds.ranges)
            transcript = TranscriptSeq.from_exons(exp.strand, cds_ref_seqs, annot.cds.ranges)

    # TODO: correct targeton region bounds based on background
    targetons = {
        s.start: s
        for s in targeton_ref_seqs
    }
    with get_db_conn() as conn:
        init_database(conn)

        if annot:

            # Load CDS annotation
            insert_exons(conn, annot.cds.ranges)

        for t in exp.targeton_configs:
            process_targeton_config(conn, codon_table, transcript, targetons, t)

        # Load background variants (targetons & exons)
        if config.bg_fp:
            # TODO: extend to support multiple contigs
            bg_ctx = get_background_context_range(targeton_ranges, annot)

            load_background_variants(conn, config.bg_fp, exp.contig, [bg_ctx])
            insert_gene_offsets(conn, bg_ctx.start, bg_ctx.end)

        if config.pam_fp:
            load_pam_variants(conn, config.pam_fp, exp.contig, targeton_ranges)

        # Load custom variants (targetons)
        if config.vcf_fp:
            load_custom_variants(conn, config.vcf_fp, exp.contig, targeton_ranges)

        # TODO: remove (DEBUG only!)
        dump_all(conn)

    # Write JSON configuration to file
    config.write(os.path.join(config.output_dir, OUTPUT_CONFIG_FILE_NAME))


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
    codon_table_fp: str | None,
    gff_fp: str | None,
    bg_fp: str | None,
    pam_fp: str | None,
    vcf_fp: str | None,

    # Output directory
    output_dir: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: str | None,
    adaptor_3: str | None,

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
