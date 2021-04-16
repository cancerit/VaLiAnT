########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
# 
# The Software remains the property of Genome Research Ltd ("GRL").
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
# 
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
# 
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
#############################

from __future__ import annotations
from collections.abc import Container, Sized
from dataclasses import dataclass
from typing import Tuple
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
        return other.strand == self.strand and super().__contains__(other)

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
