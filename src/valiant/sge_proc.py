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

import logging
import sys
from sqlite3 import Connection

from .background_variants import InvalidBackgroundVariant, RegisteredBackgroundVariant
from .codon_table import CodonTable
from .codon_table_builder import CodonTableBuilder
from .common_cli import finalise
from .constants import OUTPUT_CONFIG_FILE_NAME, OUTPUT_REF_SEQUENCES
from .contig_filter import ContigFilter
from .custom_variant import CustomVariant
from .db import get_db_conn, init_db
from .enums import SrcType
from .experiment_meta import ExperimentMeta
from .genomic_position_offsets import GenomicPositionOffsets
from .loaders.bed import BedLoader
from .loaders.experiment import ExperimentConfig, load_targetons
from .loaders.fasta import open_fasta
from .loaders.gtf import GtfLoader
from .loaders.targeton_config import TargetonConfig
from .loaders.vcf_manifest import VcfManifest
from .meta_table import MetaTable
from .oligo_generation_info import OligoGenerationInfo
from .options import Options
from .pam_variant import InvalidPamVariant, PamVariant
from .queries import insert_background_variants, insert_custom_variant_collection, insert_exons, insert_pam_protection_edits, insert_targeton_ppes, is_meta_table_empty, select_background_variant_stats, select_background_variants, select_overlapping_background_variants, select_ppe_bg_codon_overlaps, select_ppes_by_sgrna_id, clear_per_contig_tables, clear_per_targeton_tables
from .seq import Seq
from .seq_converter import apply_variants
from .sge_config import SGEConfig
from .sge_utils import fetch_sequence
from .strings.strand import Strand
from .targeton import Targeton, is_variant_frame_shifting, is_variant_nonsynonymous
from .transcript import Transcript
from .uint_range import UIntRange, UIntRangeT
from .utils import fmt_genomic_range, safe_group_by
from .variant import Variant, get_variant_duplicate_positions


def get_ctx_range(ranges: list[UIntRangeT]) -> UIntRange:
    ctx = UIntRange.span(ranges)
    if ctx.start > 1:
        # Get an extra nucleotide at the beginning
        #  (for the purposes of VCF formatting)
        ctx = ctx.offset_start(-1)
    return ctx


def load_pam_variants(conn: Connection, vcf_fp: str, contig: str, ranges: list[UIntRange]) -> None:

    # Set positional filter
    contig_ft = ContigFilter.from_ranges(contig, ranges)

    # Load variants from VCF file
    ppe_vars = PamVariant.load_vcf(vcf_fp, contig_ft)
    logging.debug("Collected %d PAM variants." % len(ppe_vars))

    # Push the variants to the database
    insert_pam_protection_edits(conn, ppe_vars)


def validate_background_variants(
    conn: Connection,
    opt: Options,
    codon_table: CodonTable,
    gpo: GenomicPositionOffsets,
    contig: str,
    strand: Strand,
    seq_ref: Seq,
    seq_bg: Seq,
    transcript: Transcript | None,
    r: UIntRange
) -> None:
    bg_vars = select_background_variants(conn, r)

    any_frame_shift = False
    any_non_syn = False

    def is_var_non_syn(v: RegisteredBackgroundVariant) -> bool:
        assert transcript
        return is_variant_nonsynonymous(
            codon_table,
            seq_ref,
            seq_bg,
            gpo,
            transcript,
            v)

    for v in bg_vars:
        if v.in_cds:
            if is_var_non_syn(v):
                any_non_syn = True
                # TODO: reassess for multiples of three
                fs = is_variant_frame_shifting(v)
                if fs:
                    any_frame_shift = True
                logging.warning(v.get_non_syn_warn_message(contig, fs))

        elif v.alt_ref_delta != 0:
            logging.warning(f"Background indel {contig}:{v} in non-coding region may impact splicing and SGE screen outcome!")

    if any_non_syn and not opt.allow_non_syn:
        raise InvalidBackgroundVariant("Invalid background: non-synonymous variants!")

    if any_frame_shift and not opt.allow_frame_shift:
        raise InvalidBackgroundVariant("Invalid background: frame-shifting variants!")


def load_background_variants(
    conn: Connection,
    vcf_fp: str,
    mask_fp: str | None,
    contig: str
) -> int:

    # Set contig filter
    # TODO: introduce a positional filter?
    contig_ft = ContigFilter.from_ranges(contig, [])

    # Load variants from VCF file
    bg_vars = CustomVariant.load_vcf(vcf_fp, contig_ft)
    logging.debug("Collected %d background variants." % len(bg_vars))

    if mask_fp:

        # Load background variant mask (ranges to exclude)
        bg_mask = BedLoader(contig).load(mask_fp)

        if bg_mask:

            # Filter out background variants based on the mask
            # TODO (performance): perform on the database
            n = len(bg_vars)
            bg_vars = [
                x for x in bg_vars
                if not any(x.pos in r for r in bg_mask)
            ]
            m = len(bg_vars)
            if m < n:
                logging.info("Discarded %d background variants due to masking." % (n - m))

    # Push the variants to the database

    insert_background_variants(conn, bg_vars)

    return len(bg_vars)


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


def lift_targeton_config(gpo: GenomicPositionOffsets, tc: TargetonConfig) -> TargetonConfig:
    bg_ref_range = gpo.ref_to_alt_range(tc.ref)
    if not bg_ref_range:
        raise ValueError("The targeton does not exist in the alternate reference!")

    bg_region_2 = gpo.ref_to_alt_range(tc.region_2)
    if not bg_region_2:
        raise ValueError("The targeton region 2 does not exist in the alternate reference!")

    return tc.alter(bg_ref_range, bg_region_2)


def get_ppe_seq(
    conn: Connection,
    seq: Seq,
    gpo: GenomicPositionOffsets | None,
    targeton: Targeton,
    targeton_alt: Targeton | None,
) -> Seq:
    if not targeton.sgrna_ids:
        return seq

    ppes = select_ppes_by_sgrna_id(conn, targeton.sgrna_ids)
    duplicate_positions = get_variant_duplicate_positions(ppes)
    if duplicate_positions:
        for pos in duplicate_positions:
            logging.warning(f"Multiple PPE's at {targeton.config.contig}:{pos}!")
        logging.critical(f"Overlapping PPE's in targeton {targeton.as_str()}!")
        sys.exit(1)

    if gpo:
        # Correct PPE coordiantes
        ppes = list(map(gpo.ref_to_alt_variant, ppes))

    try:

        # Register the [potentially] coordinate-corrected PPE's on the database
        insert_targeton_ppes(conn, ppes)

    except InvalidPamVariant:
        logging.critical(f"Multiple PPE's in the same codon in tartgeton {targeton.as_str()}!")
        sys.exit(1)

    ppes_in_range: list[Variant] = []
    for ppe in ppes:

        targeton_range = (targeton_alt or targeton).ref
        if ppe.pos in targeton_range:
            ppes_in_range.append(ppe.to_variant())
        else:
            logging.warning(
                f"Pam protection edit {ppe} outside of targeton range {targeton_range}!")

    ppes_in_range.sort(key=lambda x: x.pos)

    # Apply the PPE's to the [optionally background-corrected] sequence
    ppe_seq = apply_variants(seq, len(seq), ppes_in_range)

    return ppe_seq


def get_targeton_seq(seq: Seq, t: Targeton) -> Seq:
    return seq.subseq(t.ref, rel=False, prev_nt=True)


def proc_targeton(
    conn: Connection,
    config: SGEConfig,
    options: Options,
    gpo: GenomicPositionOffsets | None,
    ctx_seq: Seq,
    ctx_seq_bg: Seq | None,
    contig: str,
    strand: Strand,
    exp_meta: ExperimentMeta,
    exp_config: ExperimentConfig,
    codon_table: CodonTable,
    transcript_ref: Transcript | None,
    transcript_bg: Transcript | None,
    targeton_config: TargetonConfig
) -> OligoGenerationInfo:

    # Clear targeton-specific tables
    clear_per_targeton_tables.execute(conn)
    conn.commit()

    # TODO: check if this should be done once per transcript instead
    #  In theory, only upstream background variants would have an impact...?
    if gpo:
        assert ctx_seq_bg
        assert transcript_bg
        validate_background_variants(
            conn,
            options,
            codon_table,
            gpo,
            contig,
            strand,
            ctx_seq,
            ctx_seq_bg,
            transcript_bg,
            targeton_config.ref)

    targeton_ref = Targeton(targeton_config)
    targeton_alt = Targeton(lift_targeton_config(gpo, targeton_config)) if gpo else None

    ppe_seq = get_ppe_seq(
        conn,
        ctx_seq_bg or ctx_seq,
        gpo,
        targeton_ref,
        targeton_alt)

    # Test background variants against the PPE's for codon overlaps
    ppe_bg_overlap_positions = select_ppe_bg_codon_overlaps(conn)
    if ppe_bg_overlap_positions:
        for pos in ppe_bg_overlap_positions:
            logging.warning(
                f"A PAM protection edit at {contig}:{pos} overlaps a background variant in a coding region!")
        raise InvalidBackgroundVariant(
            "Invalid background: PAM protection edits overlapping background variants in a coding region!")

    # Compute the PPE mutation types
    if transcript_ref:
        ppe_mut_types = (targeton_alt or targeton_ref).get_ppe_mut_types(
            conn,
            codon_table,
            gpo,
            ctx_seq,
            ppe_seq,
            transcript_ref,
            transcript_bg or transcript_ref)
    else:
        # TODO: verify this is acceptable for non-coding!
        ppe_mut_types = []

    (targeton_alt or targeton_ref).process(
        conn, contig, codon_table, gpo, transcript_bg or transcript_ref, ppe_seq)

    targeton_name = targeton_config.name

    if not is_meta_table_empty(conn):
        t_ref = get_targeton_seq(ctx_seq, targeton_ref)
        t = targeton_alt or targeton_ref
        t_alt = get_targeton_seq(ppe_seq, t)
        if ctx_seq_bg:
            t_bg_vars = select_overlapping_background_variants(conn, t.ref)
            t_bg = get_targeton_seq(ctx_seq_bg, t)
        else:
            t_bg_vars = []
            t_bg = t_ref
        mt = MetaTable(
            SrcType.REF,
            None,
            config,
            exp_meta,
            options,
            gpo,
            contig,
            strand,
            t_ref,
            t_bg,
            t_alt,
            t_bg_vars,
            ppe_mut_types,
            transcript_ref.info if transcript_ref else None)

        # Write metadata files
        return mt.to_csv(conn, targeton_name)

    else:
        logging.warning(
            "No mutations for targeton '%s'!", targeton_name)
        return OligoGenerationInfo()


def get_gpo_ctx(
    conn: Connection,
    ctx: UIntRange
) -> tuple[GenomicPositionOffsets | None, UIntRange]:
    """
    Calculate the information required for the lift over

    Assumption: the background variants have already been loaded into the database
    """

    # TODO: also verify there are background variants in range
    bg_vars = select_background_variant_stats(conn, ctx)

    if bg_vars:

        # Extend context to liminal background variants, if necessary
        # TODO (performance): optimise (perform on the database?)
        ctx = UIntRange.span([ctx, UIntRange(
            min(x.pos for x in bg_vars),
            max(x.ref_end for x in bg_vars))
        ])
        return GenomicPositionOffsets.from_var_stats(bg_vars, ctx), ctx

    return None, ctx


def get_ctx_seq_bg(conn: Connection, gpo: GenomicPositionOffsets, ctx_seq: Seq) -> Seq:
    # TODO: fetch them only once for GPO and this...?
    bg_vars = select_background_variants(conn, ctx_seq.get_range())
    return apply_variants(ctx_seq, gpo.alt_length, bg_vars)


def proc_contig(
    conn: Connection,
    config: SGEConfig,
    exp_meta: ExperimentMeta,
    exp_config: ExperimentConfig,
    contig: str,
    codon_table_builder: CodonTableBuilder,
    targeton_configs: list[TargetonConfig],
    stats: OligoGenerationInfo
) -> None:

    # Clear the contig-specific database tables
    clear_per_contig_tables.execute(conn)

    targeton_ranges = [t.ref for t in targeton_configs]

    # Load PAM protection edits (PPE's)
    if config.pam_fp:
        load_pam_variants(conn, config.pam_fp, contig, targeton_ranges)

    # Load custom variants (targetons)
    if config.vcf_fp:
        load_custom_variants(conn, config.vcf_fp, contig, targeton_ranges)

    opt = config.get_options()
    for strand, tr_group in safe_group_by(targeton_configs, lambda x: x.strand):
        tcs: list[TargetonConfig] = list(tr_group)
        gpo: GenomicPositionOffsets | None = None
        transcript_bg: Transcript | None = None
        bg_var_count: int = 0

        # Build codon table (strand-specific)
        codon_table = codon_table_builder.build(strand)

        if config.bg_fp:

            # Load background variants
            bg_var_count = load_background_variants(
                conn, config.bg_fp, config.mask_bg_fp, contig)

        if config.gff_fp:

            # Retrieve transcript information and exon ranges
            annot = GtfLoader(contig, strand).load_gtf(config.gff_fp)
            transcript = Transcript(annot.transcript_info, annot.cds)
            assert transcript.exons

            # Compute context (preliminary)
            ctx = get_ctx_range(transcript.exons + [x.ref for x in tcs])

            if bg_var_count > 0:

                # Compute background offsets and correct the context
                gpo, ctx = get_gpo_ctx(conn, ctx)
                if gpo:

                    # Lift transcript over to background
                    transcript_bg = transcript.lift_exons(gpo)

            # TODO (future multi-transcript): clear existing exons from the database
            insert_exons(conn, strand, (transcript_bg or transcript).exons)

            with open_fasta(config.ref_fasta_fp) as fa:

                # Fetch context sequence
                ctx_seq = fetch_sequence(fa, contig, ctx)

            # Lift context sequence over to background
            ctx_seq_bg = get_ctx_seq_bg(conn, gpo, ctx_seq) if gpo else None

            for targeton_config in tcs:

                # Process targeton
                stats.update(proc_targeton(
                    conn, config, opt, gpo, ctx_seq, ctx_seq_bg, contig, strand, exp_meta, exp_config, codon_table, transcript, transcript_bg, targeton_config))

        else:
            targeton_ctxs = [
                (t, get_ctx_range([t.ref]))
                for t in tcs
            ]
            with open_fasta(config.ref_fasta_fp) as fa:
                for targeton_config, ctx in targeton_ctxs:
                    ctx = get_ctx_range([targeton_config.ref])

                    # Compute background offsets and correct the context
                    gpo, ctx = get_gpo_ctx(conn, ctx)

                    ctx_seq = fetch_sequence(fa, contig, ctx)

                    # Lift context sequence over to background
                    ctx_seq_bg = get_ctx_seq_bg(conn, gpo, ctx_seq) if gpo else None

                    stats.update(proc_targeton(
                        conn, config, opt, gpo, ctx_seq, ctx_seq_bg, contig, strand, exp_meta, exp_config, codon_table, None, None, targeton_config))


def get_ref_seq_info(seq: Seq, t: TargetonConfig) -> list[str]:
    items = [
        t.base_name,
        fmt_genomic_range(t.contig, t.ref)
    ]
    for r in t.get_all_regions():
        items.append(str(r.start))
        items.append(seq.substr(r, rel=False))
    return items


def write_ref_sequences(config: SGEConfig, targetons: list[TargetonConfig]) -> None:
    with (
        open_fasta(config.ref_fasta_fp) as fa,
        open(config.get_output_file_path(OUTPUT_REF_SEQUENCES), 'w') as fh
    ):
        for t in targetons:
            seq = fetch_sequence(fa, t.contig, t.ref)
            fh.write(','.join(get_ref_seq_info(seq, t)))
            fh.write('\n')


def run_sge(config: SGEConfig, sequences_only: bool) -> None:
    exp_meta = ExperimentMeta(
        config.species,
        config.assembly,
        config.min_length,
        config.max_length)

    # Load experiment configuration
    targeton_configs = load_targetons(config.oligo_info_fp)
    exp_config = ExperimentConfig.from_configs(targeton_configs)

    # Write reference sequences file
    write_ref_sequences(config, targeton_configs)

    if sequences_only:
        sys.exit(0)

    # Load codon table
    codon_table_builder = CodonTableBuilder.load(
        config.get_codon_table_file_path())

    # Load annotation
    if config.gff_fp and len(exp_config.targeton_configs) > 1:
        logging.warning(
            "Annotation provided for targetons on multiple contigs: " +
            "verify only one transcript is used!")

    stats = OligoGenerationInfo()

    with get_db_conn() as conn:

        # Initialise in-memory database
        init_db(conn)

        try:
            for contig, targeton_configs in exp_config.targeton_configs.items():
                proc_contig(conn, config, exp_meta, exp_config, contig, codon_table_builder, targeton_configs, stats)
        except InvalidBackgroundVariant as ex:
            logging.critical(ex.args[0])
            sys.exit(1)

    finalise(config, stats)

    # Write JSON configuration to file
    config.write(config.get_output_file_path(OUTPUT_CONFIG_FILE_NAME))
