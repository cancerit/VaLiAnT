########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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
    short_oligo_n: int
    long_oligo_n: int

    def __init__(self) -> None:
        self.short_oligo_n = 0
        self.long_oligo_n = 0

    def update(self, info: OligoGenerationInfo) -> None:
        self.short_oligo_n += info.short_oligo_n
        self.long_oligo_n += info.long_oligo_n
