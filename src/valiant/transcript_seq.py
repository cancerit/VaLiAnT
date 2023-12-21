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

from __future__ import annotations

from dataclasses import dataclass, replace

from valiant.variant import Variant
from valiant.variant_group import VariantGroup

from .cds_seq import CdsSeq
from .exon import Exon
from .seq import Seq
from .seq_collection import SeqCollection
from .strings.strand import Strand
from .uint_range import UIntRange


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

    def get_cds_seq(self, exon_index: int, r: UIntRange) -> CdsSeq:
        exon = self.get_exon(exon_index)
        ds = r.start - exon.start
        assert ds >= 0

        # TODO: remember to check whether the frame is in the GTF vs the Valiant convention
        #  (probably GTF since the loader was rewritten)
        before, after = exon.get_5p_3p_extensions(self.strand)
        return self.get_as_cds_seq(
            self.get_exon_seq_index(exon_index), r, before=before, after=after)

    def get_codon_indices(self, exon_index: int, r: UIntRange) -> list[int]:
        return self.get_exon(exon_index).get_codon_indices(r)

    def get_as_cds_seq(self, i: int, r: UIntRange, before: int = 0, after: int = 0) -> CdsSeq:
        a, b, c = self.split_substr(i, r, before=before, after=after)
        return CdsSeq(r.start, b, a, c)

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
