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
from .cds_seq import CdsSeq
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
from .variant import PatternVariant, Variant, VariantT


class InvalidTargetonRegion(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__()
        self.msg = f"Invalid targeton region: {msg}!"


def is_variant_frame_shifting(variant: Variant) -> bool:
    return variant.alt_ref_delta != 0


def get_variant_codons(transcript: Transcript, seq: Seq, variant: Variant) -> list[CdsSeq]:
    # Search all codons (supports variants spanning multiple exons)
    if variant.ref_len > 1:
        codons = transcript.get_codons_in_range(seq, variant.ref_range)
    else:
        # Handle insertion
        codon = transcript.get_codon_at(seq, variant.pos)
        codons = [codon] if codon is not None else []
    return codons


def is_variant_nonsynonymous(
    codon_table: CodonTable,
    seq_bg: Seq,
    gpo: GenomicPositionOffsets,
    variant: Variant,
    codons: list[CdsSeq]
) -> bool:
    assert codons

    if is_variant_frame_shifting(variant):
        return True

    for codon_ref in codons:
        codon_ref_range = codon_ref.ext_positions
        codon_alt_range = [gpo.ref_to_alt_position(x) for x in codon_ref_range]

        # A difference in span would imply a frame shift
        assert codon_alt_range and len(codon_alt_range) == len(codon_ref_range)

        ref = Codon(codon_ref.ext)
        alt = Codon(seq_bg.get_at(codon_ref_range))

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

    vars = get_vars_in_region(vars)
    annot_vars = get_vars_in_region(annot_vars)

    return vars, annot_vars


def _get_oligo(alt: Seq, x: VariantT, ref_start: int | None = None) -> OligoSeq[VariantT]:
    return OligoSeq.from_ref(alt, x, ref_start=ref_start)


def _get_oligos(
    variants: list[VariantT],
    get_oligo_f: Callable[[VariantT], OligoSeq[VariantT] | None]
) -> list[OligoSeq[VariantT]]:
    return [
        oligo
        for x in variants
        if (oligo := get_oligo_f(x)) is not None
    ]


@dataclass(slots=True)
class Targeton:
    config: TargetonConfig

    def as_str(self) -> str:
        return f"{self.config.contig}:{self.ref.start}-{self.ref.end}"

    @property
    def ref(self) -> UIntRange:
        return self.config.ref

    @property
    def strand(self) -> Strand:
        return self.config.strand

    @property
    def sgrna_ids(self) -> frozenset[str]:
        return self.config.sgrna_ids

    def get_mave_nt(self, variant: Variant) -> str:
        return variant.get_mave_nt(self.config.ref.start)

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
        for exon_number, ppe_ref_start, ppe_start, codon_offset in select_ppes_with_offset(conn, self.config.ref):
            rng = get_codon_range_from_offset(self.config.strand, ppe_start, codon_offset)
            assert len(rng) == 3

            codon_alt = transcript_alt.get_cds_seq(seq_alt, exon_number, rng).as_codon()
            # TODO: the original PPE start could be fetched from the database
            # TODO: verify the codon from position logic (strandedness issue?)
            codon_ref_seq = transcript_ref.get_codon_at(seq_ref, ppe_ref_start)
            assert codon_ref_seq
            codon_ref = codon_ref_seq.as_codon()
            # TODO: reconsider if this is appropriate with frame shifts

            res.append(codon_table.get_aa_change(codon_ref, codon_alt))

        return res

    def _process_custom_variants(
        self,
        conn: Connection,
        contig: str,
        gpo: GenomicPositionOffsets | None,
        targeton_seq: Seq,
    ) -> None:
        custom_vars = select_custom_variants_in_range(conn, self.config.ref)

        def get_oligo(x: VariantT) -> OligoSeq[VariantT] | None:
            if gpo:
                if gpo.ref_var_overlaps_var(x):
                    logging.warning(
                        "Custom variant %s:%s" % (contig, x) +
                        " overlaps a coordinate-shifting background variant (discarded)!")
                    return None

                # Lift custom variant positions from reference to altered
                y = gpo.ref_to_alt_variant(x)

            else:
                y = x

            return _get_oligo(targeton_seq, y, ref_start=x.pos)

        # Register the filtered custom variants
        insert_targeton_custom_variants(
            conn,
            _get_oligos(custom_vars, get_oligo),
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
        contig: str,
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

        def get_oligo(x: VariantT) -> OligoSeq[VariantT] | None:
            if gpo:
                if gpo.alt_var_overlaps_var(x):
                    logging.warning(
                        "Pattern variant %s" % self.get_mave_nt(x) +
                        " overlaps coordinate-shifting background variant (discarded)!")
                    return None

                ref_start = gpo.alt_to_ref_position(x.pos)

            else:
                ref_start = x.pos
            return _get_oligo(targeton_seq, x, ref_start=ref_start)

        def get_oligos(variants: list[VariantT]) -> list[OligoSeq[VariantT]]:
            return _get_oligos(variants, get_oligo)

        insert_pattern_variants(conn, get_oligos(pattern_variants))
        insert_annot_pattern_variants(conn, get_oligos(annot_variants))

        return pattern_variants, annot_variants

    def process(
        self,
        conn: Connection,
        contig: str,
        codon_table: CodonTable,
        gpo: GenomicPositionOffsets | None,
        transcript: Transcript | None,
        seq: Seq
    ) -> None:
        targeton_seq = seq.subseq(self.config.ref, rel=False)

        self._process_custom_variants(conn, contig, gpo, targeton_seq)
        self._process_pattern_variants(conn, contig, codon_table, gpo, transcript, seq, targeton_seq)
