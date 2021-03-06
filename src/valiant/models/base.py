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
from collections.abc import Container, Sized
from dataclasses import dataclass
from typing import Tuple, Optional
from ..utils import get_region, is_strand


@dataclass(frozen=True)
class PositionRange(Sized, Container):
    __slots__ = {'start', 'end'}

    start: int
    end: int

    def __post_init__(self) -> None:
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise TypeError("Invalid position range!")
        if self.start < 1 or self.end < self.start:
            raise ValueError("Invalid position range!")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __contains__(self, other) -> bool:
        return other.start >= self.start and other.end <= self.end

    def __eq__(self, other) -> bool:
        if not isinstance(other, PositionRange):
            raise TypeError("Unsupported operation!")

        return self.start == other.start and self.end == other.end

    # Required to allow multiprocessing to pickle the object
    def __setstate__(self, state: Tuple) -> None:
        for slot, value in state[1].items():
            object.__setattr__(self, slot, value)

    def to_tuple(self) -> Tuple[int, int]:
        return self.start, self.end

    def get_subrange_before(self, pos: int, length: Optional[int] = None) -> Optional[PositionRange]:
        if pos < 1:
            raise ValueError("Invalid position!")
        if length is None:
            if pos == 1:
                return None
            s = 1
        elif length == 0:
            return None
        else:
            if length < 1:
                raise ValueError("Invalid subsequence length!")
            s = pos - length
            if s < 1:
                raise ValueError("Invalid subsequence length: out of range!")

        return PositionRange(s, pos - 1)

    def get_subrange_after(self, pos: int, length: Optional[int] = None) -> Optional[PositionRange]:
        if pos < 1:
            raise ValueError("Invalid position!")
        if length is None:
            e = len(self)
            if pos == e:
                return None
        elif length == 0:
            return None
        else:
            if length < 1:
                raise ValueError("Invalid subsequence length!")
            e = pos + length
            if e > len(self):
                raise ValueError("Invalid subsequence length: out of range!")

        return PositionRange(pos + 1, e)


@dataclass(frozen=True)
class StrandedPositionRange(PositionRange):
    __slots__ = {'start', 'end', 'strand'}

    strand: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if not is_strand(self.strand):
            raise ValueError("Invalid strand!")

    def __eq__(self, other) -> bool:
        if not isinstance(other, StrandedPositionRange):
            raise TypeError("Unsupported operation!")
        return self.strand == other.strand and super().__eq__(other)

    def __contains__(self, other) -> bool:
        if isinstance(other, StrandedPositionRange):
            return other.strand == self.strand and super().__contains__(other)
        elif isinstance(other, PositionRange):
            return super().__contains__(other)
        else:
            raise TypeError("Unsupported operation!")

    @classmethod
    def to_plus_strand(cls, pr: PositionRange) -> StrandedPositionRange:
        return cls(pr.start, pr.end, '+')


@dataclass(frozen=True)
class GenomicPosition:
    __slots__ = {'chromosome', 'position'}

    chromosome: str
    position: int

    def __post_init__(self) -> None:
        if not self.chromosome or self.position < 1:
            raise ValueError("Invalid genomic position!")

    def __add__(self, other) -> GenomicPosition:
        if not isinstance(other, int):
            raise TypeError("Invalid position offset!")
        return GenomicPosition(self.chromosome, self.position + other)

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.position}"


@dataclass(frozen=True, init=False)
class GenomicRange(StrandedPositionRange):
    __slots__ = {'chromosome', 'start', 'end', 'strand'}

    chromosome: str

    def __init__(self, chromosome: str, start: int, end: int, strand: str) -> None:
        object.__setattr__(self, 'chromosome', chromosome)
        super().__init__(start, end, strand)

    def __post_init__(self) -> None:
        if not (
            self.chromosome
            and self.start >= 1
            and self.end >= 1
            and self.start <= self.end
            and is_strand(self.strand)
        ):
            raise ValueError("Invalid genomic range!")

    @property
    def region(self) -> str:
        return get_region(self.chromosome, self.start, self.end)

    def get_sub_range(self, start: int, end: int) -> GenomicRange:
        if start < self.start or end > self.end:
            raise ValueError("Genomic range bounds out of parent range!")
        return GenomicRange(self.chromosome, start, end, self.strand)

    def get_head(self, length: int) -> GenomicRange:
        if length <= 0:
            raise ValueError("Genomic range head length must be a strictly positive integer!")
        return self.get_sub_range(self.start, self.start + length - 1)

    def get_tail(self, length: int) -> GenomicRange:
        if length <= 0:
            raise ValueError("Genomic range tail length must be a strictly positive integer!")
        return self.get_sub_range(self.end - length + 1, self.end)

    def get_from_5_prime(self, length: int) -> GenomicRange:
        return self.get_head(length) if self.strand == '+' else self.get_tail(length)

    def get_from_3_prime(self, length: int) -> GenomicRange:
        return self.get_tail(length) if self.strand == '+' else self.get_head(length)

    def get_before_head(self, length: int) -> GenomicRange:
        return GenomicRange(
            self.chromosome,
            self.start - length,
            self.start - 1,
            self.strand)

    def get_past_tail(self, length: int) -> GenomicRange:
        return GenomicRange(
            self.chromosome,
            self.end + 1,
            self.end + length,
            self.strand)

    def get_before_5_prime(self, length: int) -> GenomicRange:
        return self.get_before_head(length) if self.strand == '+' else self.get_past_tail(length)

    def get_past_3_prime(self, length: int) -> GenomicRange:
        return self.get_past_tail(length) if self.strand == '+' else self.get_before_head(length)

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __eq__(self, other) -> bool:
        return (
            other.chromosome == self.chromosome
            and other.strand == self.strand
            and other.start == self.start
            and other.end == self.end
        )

    def __contains__(self, other) -> bool:
        return other.chromosome == self.chromosome and super().__contains__(other)

    def contains_position(self, genomic_position: GenomicPosition) -> bool:
        return (
            genomic_position.chromosome == self.chromosome
            and genomic_position.position >= self.start
            and genomic_position.position <= self.end
        )

    def as_unstranded(self) -> Tuple[str, int, int]:
        return self.chromosome, self.start, self.end

    def as_pyrange(self) -> Tuple[str, str, int, int]:
        return self.chromosome, self.strand, self.start - 1, self.end

    def as_unstranded_pyrange(self) -> Tuple[str, int, int]:
        return self.chromosome, self.start - 1, self.end

    def get_relative_subrange(self, child: GenomicRange) -> Tuple[int, int]:
        if child not in self:
            raise ValueError("Invalid subsequence!")
        start: int = child.start - self.start
        end: int = child.end - self.start + 1
        return start, end

    @property
    def pos_range(self) -> PositionRange:
        return PositionRange(self.start, self.end)


@dataclass
class TranscriptInfo:
    __slots__ = {'gene_id', 'transcript_id'}

    gene_id: str
    transcript_id: str
