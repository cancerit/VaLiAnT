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

from dataclasses import dataclass

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
        return cls(seqs, Strand(strand), {
            e.index: e
            for e in exons
        })

    def get_exon(self, exon_index: int) -> Exon:
        return self.exons[exon_index]

    def get_exon_seq_index(self, exon_index: int) -> int:
        return exon_index if self.strand.is_plus else (self.exon_count - exon_index - 1)

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
