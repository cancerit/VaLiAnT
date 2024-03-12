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

from collections.abc import Container
from dataclasses import dataclass, replace
from itertools import chain
from typing import Iterable, Sized, Generic, TypeVar


UIntRangeT = TypeVar('UIntRangeT', bound='UIntRange')


@dataclass(slots=True, frozen=True)
class UIntRange(Sized, Container):
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"Invalid range [{self.start}, {self.end}]!")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __lt__(self, other) -> bool:
        return (
            self.end < other.end if self.start == other.start else
            self.start < other.start
        )

    def __repr__(self) -> str:
        return f"[{self.start}, {self.end}]"

    def __contains__(self, x) -> bool:
        if isinstance(x, int):
            return self.start <= x <= self.end
        elif isinstance(x, UIntRange):
            return x.start in self and x.end in self
        raise TypeError("Operand type not supported!")

    def offset_start(self, offset: int):
        return replace(self, start=self.start + offset)

    def offset_end(self, offset: int):
        return replace(self, end=self.end + offset)

    def overlaps(self, other: UIntRange) -> bool:
        return (
            (other.start in self or other.end in self) or
            self in other
        )

    def intersect(self, r: UIntRange) -> UIntRange | None:
        if not self.overlaps(r):
            return None
        start = max(self.start, r.start)
        end = min(self.end, r.end)
        return UIntRange(start, end)

    @classmethod
    def from_length(cls, start: int, length: int) -> UIntRange:
        if length < 1:
            raise ValueError("Invalid range length: not strictly positive!")
        return cls(start, start + length - 1)

    @classmethod
    def from_pos(cls, pos: int) -> UIntRange:
        return UIntRange(pos, pos)

    @classmethod
    def span(cls, ranges: Iterable[UIntRangeT]) -> UIntRange:
        return cls(
            min(r.start for r in ranges),
            max(r.end for r in ranges)
        )

    @property
    def positions(self) -> list[int]:
        return list(range(self.start, self.end + 1))

    def to_slice(self, offset: int = 0) -> slice:
        return slice(self.start - offset, self.end - offset + 1)

    def to_tuple(self) -> tuple[int, int]:
        return self.start, self.end

    def offset(self, offset: int) -> UIntRange:
        return UIntRange(self.start + offset, self.end + offset)

    def get_before(self, n: int) -> UIntRange:
        assert n > 0
        return UIntRange(self.start - n, self.start - 1)

    def get_after(self, n: int) -> UIntRange:
        assert n > 0
        return UIntRange(self.end + 1, self.end + n)

    def get_until(self, pos: int) -> UIntRange:
        assert pos in self
        return replace(self, end=pos)

    def get_from(self, pos: int) -> UIntRange:
        assert pos in self
        return replace(self, start=pos)


class UIntRangeSortedList(list, Generic[UIntRangeT]):

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.sort()

    def __contains__(self, x) -> bool:
        return any(
            r.start <= x <= r.end
            for r in self
        )

    def get_before(self, i: int, r: UIntRange, before: int) -> list[int]:
        assert i >= 0 and before >= 0

        if before == 0:
            return []

        seq = self[i]
        ds = r.start - seq.start
        if before <= ds:

            # Local extension
            return list(range(r.start - before, r.start))

        else:

            # Distal & local extension
            assert i > 0
            prev_seq = self[i - 1]

            return list(chain(
                range(prev_seq.end - before + ds + 1, prev_seq.end + 1),
                range(seq.start, seq.start + ds)
            ))

    def get_after(self, i: int, r: UIntRange, after: int) -> list[int]:
        assert i >= 0 and after >= 0

        if after == 0:
            return []

        seq = self[i]
        ds = seq.end - r.end
        if after <= ds:

            # Local extension
            return list(range(r.end + 1, r.end + after + 1))

        else:

            # Distal extension
            assert i < len(self) - 1
            next_seq = self[i + 1]

            return list(chain(
                range(seq.end - ds + 1, seq.end + 1) if ds > 0 else [],
                range(next_seq.start, next_seq.start + after - ds)
            ))
