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

STRANDS = {
    '+': 'plus',
    '-': 'minus'
}


class Strand(str):
    def __init__(self, s: str) -> None:
        if s not in STRANDS:
            raise ValueError(f"Invalid strand: {s}!")
        super().__init__()

    @property
    def is_plus(self) -> bool:
        return self == '+'

    @property
    def is_minus(self) -> bool:
        return self == '-'

    @property
    def label(self) -> str:
        return STRANDS[self]
