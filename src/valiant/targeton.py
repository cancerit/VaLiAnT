########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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
from .codon_table import CodonTable
from .experiment_meta import ExperimentMeta
from .loaders.experiment import ExperimentConfig
from .loaders.targeton_config import TargetonConfig
from .meta_table import MetaTable
from .mutator import MutatorCollection
from .oligo_generation_info import OligoGenerationInfo
from .oligo_seq import OligoSeq
from .options import Options
from .pattern_variant import PatternVariant
from .queries import insert_annot_pattern_variants, insert_pattern_variants, insert_targeton_custom_variants, insert_targeton_ppes, is_meta_table_empty, select_exons_in_range, select_ppes_in_range, select_bgs_in_range, clear_per_targeton_tables, select_custom_variants_in_range
from .seq import Seq
from .sge_config import SGEConfig
from .transcript_seq import TranscriptSeq
from .uint_range import UIntRange
from .variant import RegisteredVariant, Variant, VariantT
from .variant_group import VariantGroup


GetVariantsInRangeCallable = Callable[[Connection, UIntRange], list[RegisteredVariant]]


class InvalidTargetonRegion(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__()
        self.msg = f"Invalid targeton region: {msg}!"


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


    def _fetch_variant_group(self, conn: Connection, f: GetVariantsInRangeCallable) -> VariantGroup:
        return VariantGroup.from_variants(f(conn, self.seq.get_range()))

    def apply_background_variants(self, conn: Connection) -> None:
        vars = select_bgs_in_range(conn, self.seq.get_range())
        pass

    def apply_ppes(self, conn: Connection) -> None:
        vars = select_ppes_in_range(conn, self.seq.get_range())
        pass

    def apply_custom_variants(self) -> None:
        """
        1. read all custom variants in range (background-corrected)
        2. mark those that occur in the custom regions, if any
        3. generate the oligonucleotide sequences?
        """
        pass

    def apply_pattern_variants(self) -> None:
        pass

    def apply_mutations(self) -> None:
        self.apply_custom_variants()
        self.apply_pattern_variants()

    def alter(self, conn: Connection) -> Seq:
        bg_vars = self._fetch_variant_group(conn, select_bgs_in_range)
        ppe_vars = self._fetch_variant_group(conn, select_ppes_in_range)

        bg_seq, alt_bg_vars = bg_vars.apply(self.seq, ref_check=True)
        ppe_seq, alt_ppe_vars = ppe_vars.apply(bg_seq, ref_check=False)

        # TODO: push the corrected PPE's to database
        insert_targeton_ppes(conn, alt_ppe_vars)

        return ppe_seq

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

    def _get_oligo(self, alt: Seq, x: Variant) -> OligoSeq:
        return OligoSeq.from_ref(alt, x)

    def _process_custom_variants(self, conn: Connection, opt: Options, alt: Seq) -> None:
        custom_vars = select_custom_variants_in_range(conn, self.seq.get_range())

        def get_oligo(x: Variant) -> OligoSeq:
            return self._get_oligo(alt, x)

        # Register the filtered custom variants
        insert_targeton_custom_variants(
            conn,
            list(map(get_oligo, custom_vars)),
            self.config.get_const_regions())

    def _process_pattern_variants(
        self,
        conn: Connection,
        codon_table: CodonTable,
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
            return self._get_oligo(alt, x)

        insert_pattern_variants(conn, list(map(get_oligo, pattern_variants)))
        insert_annot_pattern_variants(conn, list(map(get_oligo, annot_variants)))

        return pattern_variants, annot_variants

    def process(
        self,
        conn: Connection,
        opt: Options,
        codon_table: CodonTable,
        transcript: TranscriptSeq | None
    ) -> tuple[Seq, list[PatternVariant], list[AnnotVariant]]:

        # Truncate targeton-specific tables
        clear_per_targeton_tables.execute(conn)
        conn.commit()

        # Alter the sequence based on background variants and PPE's
        alt = self.alter(conn)
        conn.commit()

        # Process custom variants
        self._process_custom_variants(conn, opt, alt)

        # Process pattern variants
        pattern_variants, annot_variants = self._process_pattern_variants(
            conn, codon_table, transcript, alt)

        return alt, pattern_variants, annot_variants


def generate_metadata_table(
    conn: Connection,
    targeton: Targeton,
    alt: Seq,
    config: SGEConfig,
    exp: ExperimentMeta,
    exp_cfg: ExperimentConfig,
    annot: Annotation | None
) -> OligoGenerationInfo:
    targeton_name = targeton.config.name
    if not is_meta_table_empty(conn):
        meta_fn = f"{targeton_name}_meta.csv"
        meta_fp = config.get_output_file_path(meta_fn)
        options = config.get_options()
        mt = MetaTable(config, exp, options, exp_cfg, targeton.seq, alt, annot)
        return mt.to_csv(conn, meta_fp)
    else:
        logging.warning(
            "No mutations for targeton '%s'!", targeton_name)
        return OligoGenerationInfo()
