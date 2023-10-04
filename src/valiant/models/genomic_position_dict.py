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

from collections import defaultdict
from typing import Generic, TypeVar

from .base import GenomicPosition

T = TypeVar('T')


class GenomicPositionDict(Generic[T]):
    __slots__ = ['_d']

    def __init__(self) -> None:
        self._d = defaultdict(lambda: defaultdict(dict))

    def set(self, p: GenomicPosition, v: T) -> None:
        self._d[p.chromosome][p.position] = v

    def get(self, p: GenomicPosition) -> T:
        return self.get_at(p.chromosome, p.position)

    def get_at(self, chromosome: str, position: int) -> T:
        return self._d[chromosome][position]
