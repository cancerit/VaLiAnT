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

from .annotation import Exon
from .utils import get_cds_ext_3_length, get_codon_offset_complement
from .seq_collection import SeqCollection
from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .seq import Seq
from .strings.strand import Strand


@dataclass
class TranscriptSeq(SeqCollection):
    strand: Strand
    exons: dict[int, Exon]

    @classmethod
    def from_exons(cls, strand: str, seqs: list[Seq], exons: list[Exon]) -> TranscriptSeq:
        return cls(seqs, Strand(strand), {
            e.index: e
            for e in exons
        })

    def get_exon(self, exon_index: int) -> Exon:
        return self.exons[exon_index]

    def get_codon_offset(self, exon_index: int, pos: int) -> int:
        exon = self.get_exon(exon_index)
        return pos - exon.start + exon.frame

    def get_ext(self, exon_index: int, r: UIntRange) -> DnaStr:
        exon = self.get_exon(exon_index)
        ds = r.start - exon.start
        assert ds >= 0

        cds_prefix_len = get_codon_offset_complement(exon.frame)
        cds_suffix_len = get_cds_ext_3_length(cds_prefix_len, len(exon))

        before, after = (
            (cds_prefix_len, cds_suffix_len) if self.strand.is_plus else
            (cds_suffix_len, cds_prefix_len)
        )
        return self.substr(exon_index, r, before=before, after=after)
