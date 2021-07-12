########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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
from typing import Dict, Optional, Iterable, FrozenSet, List, Tuple
from ..loaders.fasta import load_from_multi_fasta
from ..loaders.tsv import load_tsv
from .base import PositionRange, StrandedPositionRange
from .cdna import CDNA, AnnotatedCDNA
from .sequences import Sequence
from .sequence_info import SequenceInfo


CDNA_ANNOT_FIELDS = [
    'seq_id',
    'gene_id',
    'transcript_id',
    'cds_start',
    'cds_end'
]


def load_seqs(fp: str, ids: Iterable[str]) -> Dict[str, Sequence]:
    return {
        seq_id: Sequence(seq)
        for seq_id, seq in load_from_multi_fasta(fp, ids).items()
    }


def load_annot(fp: str, ids: Iterable[str]) -> Dict[str, Tuple[SequenceInfo, Optional[StrandedPositionRange]]]:
    def parse_cds_range(cds_start: str, cds_end: str) -> Optional[StrandedPositionRange]:
        if bool(cds_start) != bool(cds_end):
            raise ValueError("Invalid CDS range!")
        if not cds_start:
            return None
        return StrandedPositionRange(
            int(cds_start), int(cds_end), '+')

    return {
        seq_id: (
            SequenceInfo(gene_id or None, transcript_id or None),
            parse_cds_range(cds_start, cds_end)
        )
        for seq_id, gene_id, transcript_id, cds_start, cds_end in load_tsv(
            fp, CDNA_ANNOT_FIELDS)
        if seq_id in ids
    }


# TODO: currently supporting only all-or-none cDNA annotation...?
@dataclass
class CDNASequenceRepository:
    __slots__ = {'_sequences'}

    _sequences: Dict[str, CDNA]

    @classmethod
    def load(cls, ids: FrozenSet[str], fasta_fp: str, annot_fp: Optional[str] = None) -> CDNASequenceRepository:
        def get_cdna(seq: Sequence, info: SequenceInfo, cds: Optional[StrandedPositionRange]) -> CDNA:
            return AnnotatedCDNA(seq, info, cds) if cds else CDNA(seq, info)

        seq_id_seqs: Dict[str, Sequence] = load_seqs(fasta_fp, ids)
        if annot_fp:
            seq_id_cds = load_annot(annot_fp, ids)
            return cls({
                seq_id: get_cdna(seq_id_seqs[seq_id], *seq_id_cds[seq_id])
                for seq_id in ids
            })
        else:
            return cls({
                seq_id: CDNA(seq_id_seqs[seq_id], SequenceInfo.empty())
                for seq_id in ids
            })

    def get(self, seq_id: str) -> Optional[CDNA]:
        return self._sequences.get(seq_id, None)

    def get_subsequence(self, seq_id: str, pr: PositionRange) -> Optional[Sequence]:
        if seq_id not in self._sequences:
            return None

        cdna = self._sequences[seq_id]
        return cdna.get_subsequence(pr)
