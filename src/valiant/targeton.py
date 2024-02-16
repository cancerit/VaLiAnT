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
from typing import Iterable

from .annot_variant import AnnotVariant
from .background_variants import RegisteredBackgroundVariant
from .codon_table import CodonTable
from .enums import MutationType
from .exon import get_codon_range_from_offset
from .genomic_position_offsets import GenomicPositionOffsets
from .loaders.targeton_config import TargetonConfig
from .mutator import MutatorCollection
from .oligo_seq import OligoSeq
from .queries import insert_annot_pattern_variants, insert_pattern_variants, insert_targeton_custom_variants, select_exons_in_range, select_ppes_with_offset, select_custom_variants_in_range
from .seq import Seq
from .strings.codon import Codon
from .strings.strand import Strand
from .transcript import Transcript
from .uint_range import UIntRange
from .variant import PatternVariant, RegisteredVariant, Variant, VariantT


class InvalidTargetonRegion(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__()
        self.msg = f"Invalid targeton region: {msg}!"


def is_variant_frame_shifting(variant: RegisteredVariant) -> bool:
    return variant.alt_ref_delta != 0


def is_variant_nonsynonymous(
    strand: Strand,
    codon_table: CodonTable,
    seq: Seq,
    seq_bg: Seq,
    gpo: GenomicPositionOffsets,
    transcript: Transcript,
    variant: RegisteredBackgroundVariant
) -> bool:
    assert variant.in_cds

    if is_variant_frame_shifting(variant):
        return True

    # TODO: drop these once the new approach is validated...
    if bool(variant.start_exon_index) != bool(variant.end_exon_index):
        logging.warning("[LIMITATION] Background variant spanning coding and non-coding regions, affecting more than a single codon: may or may not be synonymous (assuming it is by default)!")
        return False

    if variant.start_exon_index != variant.end_exon_index:
        logging.warning("[LIMITATION] Background variant spanning more than one exon: may or may not be synonymous (assuming it is by default)!")
        return False

    assert variant.start_exon_index
    assert variant.start_codon_index
    assert variant.end_codon_index

    # TODO: use GPO to convert the codon coordinates to fetch ALT's

    if (
        variant.start_exon_index == variant.end_exon_index and
        variant.start_codon_index == variant.end_codon_index
    ):

        # Fall back to single codon index
        codons = [transcript.get_codon(
            seq, variant.start_exon_index, variant.start_codon_index)]

    else:

        # Search all codons (supports variants spanning multiple exons)
        codons = transcript.get_codons_in_range(seq, variant.ref_range)

    for codon_ref in codons:
        codon_ref_range = codon_ref.get_range()
        codon_alt_range = gpo.ref_to_alt_range(codon_ref.get_range())

        # A difference in span would imply a frame shift
        assert codon_alt_range and len(codon_alt_range) == len(codon_ref_range)

        # Reference transcript applied to background sequence (?)
        codon_alt = transcript.get_codon_at(seq_bg, codon_alt_range.start)
        assert codon_alt

        ref = Codon(codon_ref.ext)
        alt = Codon(codon_alt.ext)
        if not codon_table.is_syn(ref, alt):
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
    transcript: Transcript | None,
    seq: Seq,
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
        r_seq = transcript.get_cds_seq(seq, exon_id, region)
    else:
        r_seq = seq.subseq(region, rel=False)
    assert r_seq.start == region.start

    vars, annot_vars = mc.get_variants(codon_table, r_seq)

    # if is_cds:
    vars = get_vars_in_region(vars)
    annot_vars = get_vars_in_region(annot_vars)

    return vars, annot_vars


def _get_oligo(alt: Seq, x: Variant, ref_start: int | None = None) -> OligoSeq:
    return OligoSeq.from_ref(alt, x, ref_start=ref_start)


@dataclass(slots=True)
class Targeton:
    config: TargetonConfig

    @property
    def ref(self) -> UIntRange:
        return self.config.ref

    @property
    def strand(self) -> Strand:
        return self.config.strand

    @property
    def sgrna_ids(self) -> frozenset[str]:
        return self.config.sgrna_ids

    def get_ppe_mut_types(
        self,
        conn: Connection,
        codon_table: CodonTable,
        gpo: GenomicPositionOffsets | None,
        seq_ref: Seq,
        seq_alt: Seq,
        transcript_ref: Transcript,
        transcript_alt: Transcript
    ) -> list[MutationType]:
        if not self.config.sgrna_ids:
            return []

        res = []
        for exon_number, ppe_start, codon_offset in select_ppes_with_offset(conn, self.config.ref):
            rng = get_codon_range_from_offset(self.config.strand, ppe_start, codon_offset)
            assert len(rng) == 3

            codon_alt = transcript_alt.get_cds_seq(seq_alt, exon_number, rng).as_codon()
            # TODO: the original PPE start could be fetched from the database
            ppe_ref_start = gpo.alt_to_ref_position(ppe_start) if gpo else ppe_start
            assert ppe_ref_start is not None
            # TODO: verify the codon from position logic (strandedness issue?)
            # codon_ref_seq = transcript_ref.get_codon_at(seq_ref, ppe_ref_start)
            # assert codon_ref_seq
            # codon_ref = codon_ref_seq.as_codon()
            # TODO: reconsider if this is appropriate with frame shifts
            codon_ref = transcript_ref.get_cds_seq(seq_ref, exon_number, rng).as_codon()

            res.append(codon_table.get_aa_change(codon_ref, codon_alt))

        return res

    def _process_custom_variants(
        self,
        conn: Connection,
        gpo: GenomicPositionOffsets | None,
        targeton_seq: Seq,
        ref_start: int
    ) -> None:
        custom_vars = select_custom_variants_in_range(conn, self.config.ref)

        def get_oligo(x: Variant) -> OligoSeq:
            # Lift custom variant positions from reference to altered
            # TODO: filter out on clash!
            return _get_oligo(
                targeton_seq, gpo.ref_to_alt_variant(x) if gpo else x,
                ref_start=ref_start)

        # Register the filtered custom variants
        insert_targeton_custom_variants(
            conn,
            list(map(get_oligo, custom_vars)),
            self.config.get_const_regions())

    def _process_region(
        self,
        conn: Connection,
        codon_table: CodonTable,
        transcript: Transcript | None,
        seq: Seq,
        r: UIntRange,
        mc: MutatorCollection
    ) -> tuple:
        pattern_variants, annot_variants = get_pattern_variants_from_region(
            conn, codon_table, transcript, seq, r, mc)

        return pattern_variants, annot_variants

    @property
    def mutable_regions(self) -> list[tuple[UIntRange, MutatorCollection]]:
        return [
            (r, MutatorCollection.from_configs(m))
            for r, m in self.config.get_mutable_regions()
        ]

    def _process_pattern_variants(
        self,
        conn: Connection,
        codon_table: CodonTable,
        gpo: GenomicPositionOffsets | None,
        transcript: Transcript | None,
        seq: Seq,
        targeton_seq: Seq
    ) -> tuple[list[PatternVariant], list[AnnotVariant]]:
        pattern_variants: list[PatternVariant] = []
        annot_variants: list[AnnotVariant] = []

        process_region_f = partial(
            self._process_region, conn, codon_table, transcript, seq)

        for r, mc in self.mutable_regions:
            vars, annot_vars = process_region_f(r, mc)

            pattern_variants.extend(vars)
            annot_variants.extend(annot_vars)

        def get_oligo(x: Variant) -> OligoSeq:
            # TODO: check that it is appropriate to use the variant position as ref start!
            ref_start = gpo.alt_to_ref_position(x.pos) if gpo else x.pos
            return _get_oligo(targeton_seq, x, ref_start=ref_start)

        insert_pattern_variants(conn, list(map(get_oligo, pattern_variants)))
        insert_annot_pattern_variants(conn, list(map(get_oligo, annot_variants)))

        return pattern_variants, annot_variants

    def process(
        self,
        conn: Connection,
        codon_table: CodonTable,
        gpo: GenomicPositionOffsets | None,
        transcript: Transcript | None,
        seq: Seq
    ) -> None:
        targeton_seq = seq.subseq(self.config.ref, rel=False)
        # Verify it is appropriate to set it for custom variants!
        ref_start = seq.start

        self._process_custom_variants(conn, gpo, targeton_seq, ref_start)
        self._process_pattern_variants(conn, codon_table, gpo, transcript, seq, targeton_seq)
