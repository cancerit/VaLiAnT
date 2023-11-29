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

from dataclasses import dataclass
from sqlite3 import Connection

from .annot_variant import AnnotVariant
from .codon_table import CodonTable
from .loaders.experiment import TargetonConfig
from .mutator import MutatorCollection
from .pattern_variant import PatternVariant
from .queries import insert_annot_pattern_variants, insert_pattern_variants, select_exons_in_range
from .seq import Seq
from .transcript_seq import TranscriptSeq
from .uint_range import UIntRange


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
    # TODO: SNV's should be annotated for amino acid changes
    # TODO: assign partial codons when interacting with PPE's
    # TODO: assign sgRNA ID's based on PPE's

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


@dataclass(slots=True)
class Targeton:
    seq: Seq
    config: TargetonConfig

    def process(
        self,
        conn: Connection,
        codon_table: CodonTable,
        transcript: TranscriptSeq | None
    ) -> tuple[list[PatternVariant], list[AnnotVariant]]:
        pattern_variants: list[PatternVariant] = []
        annot_variants: list[AnnotVariant] = []
        t = self.config

        # TODO: apply pattern variants to R1 and R3 as well
        for i, r in enumerate([None, t.region_2, None]):
            if r is not None:
                mutators = t.mutators[i]
                if mutators:
                    mc = MutatorCollection.from_configs(mutators)
                    vars, annot_vars = get_pattern_variants_from_region(
                        conn, codon_table, transcript, self.seq, r, mc)
                    pattern_variants.extend(vars)
                    annot_variants.extend(annot_vars)

        return pattern_variants, annot_variants
