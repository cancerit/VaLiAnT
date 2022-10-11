########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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
from typing import Optional, Sized, Tuple
from .base import PositionRange, StrandedPositionRange
from .sequences import Sequence
from .sequence_info import SequenceInfo


@dataclass(frozen=True)
class CDNA(Sized):
    __slots__ = {'seq', 'seq_info'}

    seq: Sequence
    seq_info: SequenceInfo

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def range(self) -> StrandedPositionRange:
        return StrandedPositionRange(1, len(self), '+')

    def get_subsequence(self, pr: PositionRange) -> Sequence:
        return self.seq.get_rel_subsequence(pr)

    def get_subsequence_string(self, pr: PositionRange) -> str:
        return self.get_subsequence(pr).sequence

    def get_subsequence_string_before(self, pos: int, length: Optional[int] = None) -> str:
        r = self.range.get_subrange_before(pos, length=length)
        return self.get_subsequence_string(r) if r is not None else ''

    def get_subsequence_string_after(self, pos: int, length: Optional[int] = None) -> str:
        r = self.range.get_subrange_after(pos, length=length)
        return self.get_subsequence_string(r) if r is not None else ''


@dataclass(frozen=True)
class AnnotatedCDNA(CDNA):
    __slots__ = {'seq', 'seq_info', 'cds_range'}

    cds_range: StrandedPositionRange

    def _is_cds_overlap_valid(self, pr: PositionRange) -> bool:
        cds_start: int = self.cds_range.start
        if pr.start < cds_start:
            if pr.end >= cds_start:
                return False
        elif pr.end > self.cds_range.end:
            return False
        return True

    def get_extended_subsequence(self, pr: PositionRange) -> Tuple[Sequence, str, str]:
        if not self._is_cds_overlap_valid(pr):
            raise ValueError("Position range spans CDS and non-CDS regions!")

        # Get subsequence
        seq: Sequence = self.get_subsequence(pr)

        # Get 5' CDS extension
        delta_start: int = pr.start - self.cds_range.start
        ext_5p_len: int = delta_start % 3
        ext_5p: str = self.get_subsequence_string_before(
            pr.start, ext_5p_len) if ext_5p_len != 0 else ''

        # Get 3' CDS extension
        ext_3p_len: int = (3 - (len(pr) + ext_5p_len) % 3) % 3
        ext_3p: str = self.get_subsequence_string_after(
            pr.end, ext_3p_len) if ext_3p_len != 0 else ''

        return seq, ext_5p, ext_3p

    def __post_init__(self) -> None:
        n: int = len(self)
        if len(self.cds_range) % 3:
            raise ValueError("Invalid CDS length!")
        if self.cds_range.start > n or self.cds_range.end > n:
            raise ValueError("CDS out of range!")
