########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023, 2024 Genome Research Ltd
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


@dataclass(slots=True, init=False)
class OligoGenerationInfo:
    too_short: int
    in_range: int
    too_long: int

    @property
    def short_oligo_n(self) -> int:
        return self.too_short

    @property
    def long_oligo_n(self) -> int:
        return self.too_long

    @property
    def out_of_range_n(self) -> int:
        return self.short_oligo_n + self.long_oligo_n

    def __init__(self) -> None:
        self.too_short = 0
        self.in_range = 0
        self.too_long = 0

    def update(self, info: OligoGenerationInfo) -> None:
        self.too_short += info.too_short
        self.in_range += info.in_range
        self.too_long += info.too_long
