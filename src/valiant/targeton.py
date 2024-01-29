########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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
from dataclasses import dataclass
from functools import partial
from sqlite3 import Connection
from typing import Callable, Iterable

from .annot_variant import AnnotVariant
from .annotation import Annotation
from .background_variants import InvalidBackgroundVariant, RegisteredBackgroundVariant
from .codon_table import CodonTable
from .enums import MutationType
from .exon import get_codon_range, get_codon_range_from_offset
from .experiment_meta import ExperimentMeta
from .genomic_position_offsets import GenomicPositionOffsets
from .loaders.experiment import ExperimentConfig
from .loaders.targeton_config import TargetonConfig
from .meta_table import MetaTable
from .mutator import MutatorCollection
from .oligo_generation_info import OligoGenerationInfo
from .oligo_seq import OligoSeq
from .options import Options
from .pattern_variant import PatternVariant
from .queries import insert_annot_pattern_variants, insert_pattern_variants, insert_targeton_custom_variants, insert_targeton_ppes, is_meta_table_empty, select_background_variants, select_exons_in_range, clear_per_targeton_tables, select_custom_variants_in_range, select_ppe_bg_codon_overlaps, select_ppes_with_offset, sql_select_ppes_in_range
from .seq import Seq
from .sge_config import SGEConfig
from .sql_gen import SqlQuery, sql_and, sql_eq_or_in_str_list
from .strings.codon import Codon
from .strings.strand import Strand
from .transcript_seq import TranscriptSeq
from .uint_range import UIntRange
from .variant import RegisteredVariant, Variant, VariantT
from .variant_group import VariantGroup
from .variant_select import VariantSelectStart


GetVariantsInRangeCallable = Callable[[Connection, UIntRange], list[RegisteredVariant]]


class InvalidTargetonRegion(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__()
        self.msg = f"Invalid targeton region: {msg}!"


def is_variant_frame_shifting(variant: RegisteredVariant) -> bool:
    return variant.alt_ref_delta != 0


def is_variant_nonsynonymous(
    strand: Strand,
    codon_table: CodonTable,
    transcript: TranscriptSeq,
    variant: RegisteredBackgroundVariant
) -> bool:
    assert variant.in_cds

    if is_variant_frame_shifting(variant):
        return True

    if bool(variant.start_exon_index) != bool(variant.end_exon_index):
        logging.warning("[LIMITATION] Background variant spanning coding and non-coding regions, affecting more than a single codon: may or may not be synonymous (assuming it is by default)!")
        return False

    if variant.start_exon_index != variant.end_exon_index:
        logging.warning("[LIMITATION] Background variant spanning more than one exon: may or may not be synonymous (assuming it is by default)!")
        return False

    assert variant.start_exon_index
    assert variant.start_codon_index
    assert variant.end_codon_index

    exon_index = variant.start_exon_index
    alt = transcript.alter_exon(exon_index, variant)
    exon = transcript.get_exon(exon_index)
    for codon_index in range(variant.start_codon_index, variant.end_codon_index + 1):
        ref_codon = Codon(transcript.get_codon(exon_index, codon_index).s)
        origin = exon.get_first_codon_start(strand)
        alt_codon = Codon(alt.ext_substr(
            get_codon_range(strand, origin, codon_index), rel=False))
        if codon_table.get_aa_change(ref_codon, alt_codon) == MutationType.SYNONYMOUS:
            return True

    return False


def get_targeton_region_exon_id(conn: Connection, r: UIntRange) -> int | None:
    exons = select_exons_in_range(conn, r.start, r.end)
    match len(exons):
        case 0:
            return None
        case 1:
            exon_id, exon_range = exons[0]
            if r not in exon_range:
                raise InvalidTargetonRegion("contains both exonic and intronic sequences")
            return exon_id
        case _:
            raise InvalidTargetonRegion("overlaps multiple exons")


def get_pattern_variants_from_region(
    conn: Connection,
    codon_table: CodonTable,
    transcript: TranscriptSeq | None,
    targeton: Seq,
    region: UIntRange,
    mc: MutatorCollection
) -> tuple[list[PatternVariant], list[AnnotVariant]]:

    def get_vars_in_region(a: Iterable[VariantT]) -> list[VariantT]:
        return [x for x in a if x.pos in region and x.ref_end in region]

    # Get overlapping exon ID
    exon_id = get_targeton_region_exon_id(conn, region)
    is_cds = exon_id is not None

    if is_cds:
        assert transcript is not None
        r_seq = transcript.get_cds_seq(exon_id, region)
    else:
        r_seq = targeton.subseq(region, rel=False)

    vars, annot_vars = mc.get_variants(codon_table, r_seq)

    if is_cds:
        vars = get_vars_in_region(vars)
        annot_vars = get_vars_in_region(annot_vars)

    return vars, annot_vars


@dataclass(slots=True)
class Targeton:
    seq: Seq
    config: TargetonConfig

    @property
    def select_ppes_in_range(self) -> GetVariantsInRangeCallable:

        # Build query (add sgrna ID filter to positional filter)
        query = SqlQuery(sql_and([
            sql_select_ppes_in_range,
            sql_eq_or_in_str_list('sgrna_id', list(self.config.sgrna_ids))
        ]))

        return VariantSelectStart(query).select_in_range

    def fetch_background_variants(self, conn: Connection, contig: str, codon_table: CodonTable, transcript: TranscriptSeq | None, opt: Options) -> list[RegisteredBackgroundVariant]:
        bg_vars = select_background_variants(conn, self.seq.get_range())

        any_frame_shift = False
        any_non_syn = False

        for v in bg_vars:
            if v.in_cds:
                assert transcript
                if is_variant_nonsynonymous(self.config.strand, codon_table, transcript, v):
                    any_non_syn = True
                    # TODO: reassess for multiples of three
                    fs = is_variant_frame_shifting(v)
                    if fs:
                        any_frame_shift = True
                    logging.warning(v.get_non_syn_warn_message(contig, fs))

            elif v.alt_ref_delta != 0:
                logging.warning("Background indel in non-coding region may impact splicing and SGE screen outcome!")

        if any_non_syn and not opt.allow_non_syn:
            raise InvalidBackgroundVariant("Invalid background: non-synonymous variants!")

        if any_frame_shift and not opt.allow_frame_shift:
            raise InvalidBackgroundVariant("Invalid background: frame-shifting variants!")

        return bg_vars

    def alter(self, conn: Connection, contig: str, bg_vars: list[RegisteredBackgroundVariant]) -> tuple[Seq, GenomicPositionOffsets | None]:

        # Apply background variants
        if bg_vars:
            bg_seq, _ = VariantGroup.from_variants(bg_vars).apply(self.seq, ref_check=True)
            gpo = GenomicPositionOffsets.from_variants(self.seq.start, len(self.seq), bg_vars)
        else:
            bg_seq = self.seq
            gpo = None

        # Compute the genomic position offsets
        # TODO: because the background variants may be upstream, the targeton start position itself may need shifting!
        #   (This may result in negative offsets, unlikely to function properly...)

        if self.config.sgrna_ids:

            # Fetch PPE's
            ppe_vars = self.select_ppes_in_range(conn, self.seq.get_range())

            # Lift PPE's positions from reference to altered
            alt_ppe_vars = VariantGroup.from_variants(
                list(map(gpo.ref_to_alt_variant, ppe_vars)) if gpo else
                ppe_vars)

            ppe_seq, _ = alt_ppe_vars.apply(bg_seq, ref_check=False)

            # Register the coordinate-corrected PPE's on the database
            insert_targeton_ppes(conn, alt_ppe_vars.variants)

            # Test background variants against the PPE's for codon overlaps
            ppe_bg_overlap_positions = select_ppe_bg_codon_overlaps(conn)
            if ppe_bg_overlap_positions:
                for pos in ppe_bg_overlap_positions:
                    logging.warning(f"A PAM protection edit at {contig}:{pos} overlaps a background variant in a coding region!")
                raise InvalidBackgroundVariant("Invalid background: PAM protection edits overlapping background variants in a coding region!")

        else:

            # No PPE's
            ppe_seq = bg_seq

        return ppe_seq, gpo

    def _process_region(
        self,
        conn: Connection,
        codon_table: CodonTable,
        transcript: TranscriptSeq | None,
        targeton_seq: Seq,
        r: UIntRange,
        mc: MutatorCollection
    ) -> tuple:
        pattern_variants, annot_variants = get_pattern_variants_from_region(
            conn, codon_table, transcript, targeton_seq, r, mc)

        return pattern_variants, annot_variants

    @property
    def mutable_regions(self) -> list[tuple[UIntRange, MutatorCollection]]:
        return [
            (r, MutatorCollection.from_configs(m))
            for r, m in self.config.get_mutable_regions()
        ]

    def _get_oligo(self, alt: Seq, x: Variant, ref_start: int | None = None) -> OligoSeq:
        return OligoSeq.from_ref(alt, x, ref_start=ref_start)

    def _process_custom_variants(self, conn: Connection, gpo: GenomicPositionOffsets | None, alt: Seq) -> None:
        custom_vars = select_custom_variants_in_range(conn, self.seq.get_range())

        def get_oligo(x: Variant) -> OligoSeq:
            # Lift custom variant positions from reference to altered
            return self._get_oligo(alt, gpo.ref_to_alt_variant(x) if gpo else x)

        # Register the filtered custom variants
        insert_targeton_custom_variants(
            conn,
            list(map(get_oligo, custom_vars)),
            self.config.get_const_regions())

    def _process_pattern_variants(
        self,
        conn: Connection,
        codon_table: CodonTable,
        gpo: GenomicPositionOffsets | None,
        transcript: TranscriptSeq | None,
        alt: Seq
    ) -> tuple[list[PatternVariant], list[AnnotVariant]]:
        pattern_variants: list[PatternVariant] = []
        annot_variants: list[AnnotVariant] = []

        process_region_f = partial(
            self._process_region, conn, codon_table, transcript, alt)

        for r, mc in self.mutable_regions:
            vars, annot_vars = process_region_f(r, mc)

            pattern_variants.extend(vars)
            annot_variants.extend(annot_vars)

        def get_oligo(x: Variant) -> OligoSeq:
            ref_start = gpo.alt_to_ref_position(x.pos) if gpo else x.pos
            return self._get_oligo(alt, x, ref_start=ref_start)

        insert_pattern_variants(conn, list(map(get_oligo, pattern_variants)))
        insert_annot_pattern_variants(conn, list(map(get_oligo, annot_variants)))

        return pattern_variants, annot_variants

    def _process_ppe_codon(
        self,
        codon_table: CodonTable,
        transcript_bg: TranscriptSeq,
        transcript_ppe: TranscriptSeq,
        exon_index: int,
        ppe_start: int,
        codon_offset: int
    ) -> MutationType:
        rng = get_codon_range_from_offset(self.config.strand, ppe_start, codon_offset)
        assert len(rng) == 3

        codon_ref = transcript_bg.get_cds_seq(exon_index, rng).as_codon()
        codon_alt = transcript_ppe.get_cds_seq(exon_index, rng).as_codon()

        return codon_table.get_aa_change(codon_ref, codon_alt)

    def process(
        self,
        conn: Connection,
        gpo: GenomicPositionOffsets | None,
        codon_table: CodonTable,
        alt: Seq,
        transcript_bg: TranscriptSeq | None,
        transcript_ppe: TranscriptSeq | None
    ) -> tuple[list[MutationType], list[PatternVariant], list[AnnotVariant]]:

        assert bool(transcript_bg) == bool(transcript_ppe)

        # Categorise PPE's (synonymous, missense, or nonsense)
        ppe_mut_types = [
            self._process_ppe_codon(
                codon_table,
                transcript_bg,
                transcript_ppe,
                exon_index,
                ppe_start,
                codon_offset)
            for exon_index, ppe_start, codon_offset
            in select_ppes_with_offset(conn, self.seq.get_range())
        ] if self.config.sgrna_ids and transcript_bg and transcript_ppe else []

        # Process custom variants
        self._process_custom_variants(conn, gpo, alt)

        # Process pattern variants
        pattern_variants, annot_variants = self._process_pattern_variants(
            conn, codon_table, gpo, transcript_ppe, alt)

        return ppe_mut_types, pattern_variants, annot_variants


def generate_metadata_table(
    conn: Connection,
    targeton: Targeton,
    alt: Seq,
    ppe_mut_types: list[MutationType],
    gpo: GenomicPositionOffsets | None,
    config: SGEConfig,
    exp: ExperimentMeta,
    exp_cfg: ExperimentConfig,
    annot: Annotation | None
) -> OligoGenerationInfo:
    targeton_name = targeton.config.name
    if not is_meta_table_empty(conn):
        options = config.get_options()
        mt = MetaTable(config, exp, options, gpo, exp_cfg, targeton.seq, alt, ppe_mut_types, annot)
        return mt.to_csv(conn, targeton_name)
    else:
        logging.warning(
            "No mutations for targeton '%s'!", targeton_name)
        return OligoGenerationInfo()
