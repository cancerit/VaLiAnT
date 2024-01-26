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

from __future__ import annotations

from dataclasses import dataclass, replace

from .cds_seq import CdsSeq
from .codon_table import CodonTable
from .exon import Exon
from .oligo_seq import alter_seq
from .seq import Seq
from .seq_collection import SeqCollection
from .strings.codon import Codon
from .strings.strand import Strand
from .uint_range import UIntRange
from .utils import get_cds_ext_3_length
from .variant import Variant
from .variant_group import VariantGroup


def get_range_cds_exts(strand: Strand, exon: Exon, r: UIntRange) -> tuple[int, int]:
    delta_5p = (
        (r.start - exon.start) if strand.is_plus else
        (exon.end - r.end)
    )
    assert delta_5p >= 0

    cds_prefix_length = (exon.cds_prefix_length + (delta_5p % 3)) % 3
    cds_suffix_length = get_cds_ext_3_length(cds_prefix_length, len(r))

    return (
        (cds_prefix_length, cds_suffix_length) if strand.is_plus else
        (cds_suffix_length, cds_prefix_length)
    )


@dataclass(slots=True)
class TranscriptSeq(SeqCollection):
    strand: Strand
    exons: dict[int, Exon]

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @classmethod
    def from_exons(cls, strand: str, seqs: list[Seq], exons: list[Exon]) -> TranscriptSeq:
        assert len(exons) == len(seqs)
        return cls(seqs, Strand(strand), {
            e.index: e
            for e in exons
        })

    def begins_with_start_codon(self, codon_table: CodonTable) -> bool:
        first_exon = self.seqs[self.get_exon_seq_index(0)].s
        return (
            first_exon.startswith(codon_table.start_codon) if self.strand.is_plus else
            first_exon.endswith(codon_table.start_codon)
        )

    def get_exon(self, exon_index: int) -> Exon:
        return self.exons[exon_index]

    def _get_compl_index(self, i: int) -> int:
        return len(self.seqs) - i - 1

    def get_exon_seq_index(self, exon_index: int) -> int:
        return exon_index if self.strand.is_plus else self._get_compl_index(exon_index)

    def get_exon_index_from_seq_index(self, seq_index: int) -> int:
        return seq_index if self.strand.is_plus else self._get_compl_index(seq_index)

    def get_codon_offset(self, exon_index: int, pos: int) -> int:
        return self.get_exon(exon_index).get_codon_offset(self.strand, pos)

    def get_codon(self, exon_index: int, codon_index: int) -> Seq:
        exon = self.get_exon(exon_index)
        codon_range = exon.get_codon(self.strand, codon_index)
        cds_seq = self.get_cds_seq(exon_index, codon_range)
        assert not cds_seq.cds_prefix and not cds_seq.cds_suffix
        return Seq(cds_seq.start, Codon(cds_seq.s))

    def get_cds_seq(self, exon_index: int, r: UIntRange) -> CdsSeq:
        exon = self.get_exon(exon_index)

        if r not in exon:
            raise ValueError("CDS range can't span more than one exon!")

        before, after = get_range_cds_exts(self.strand, exon, r)

        return self.get_as_cds_seq(
            self.get_exon_seq_index(exon_index), r, before=before, after=after)

    def get_codon_indices(self, exon_index: int, r: UIntRange) -> list[int]:
        return self.get_exon(exon_index).get_codon_indices(r)

    def get_as_cds_seq(self, i: int, r: UIntRange, before: int = 0, after: int = 0) -> CdsSeq:
        a, b, c = self.split_substr(i, r, before=before, after=after)
        return CdsSeq(r.start, b, cds_prefix=a, cds_suffix=c)

    def alter(self, exon_ppes: dict[int, VariantGroup[Variant]]) -> TranscriptSeq:
        # TODO: apply background variants as well

        def alter_exon(seq_index: int) -> Seq:
            exon_index = self.get_exon_index_from_seq_index(seq_index)
            seq = self.seqs[seq_index]

            return (
                exon_ppes[exon_index].apply_no_offset(
                    self.seqs[seq_index]) if exon_index in exon_ppes else
                seq
            )

        alt_seqs = [
            alter_exon(i)
            for i in range(len(self.seqs))
        ]

        # TODO: alter the exon ranges as well after applying background variants
        return replace(self, seqs=alt_seqs)

    def get_exon_seq(self, exon_index: int) -> CdsSeq:
        return self.get_cds_seq(exon_index, self.get_exon(exon_index))

    def alter_exon(self, exon_index: int, variant: Variant) -> CdsSeq:
        ref = self.get_exon_seq(exon_index)
        alt = alter_seq(ref, variant)
        return replace(ref, s=alt.s)
