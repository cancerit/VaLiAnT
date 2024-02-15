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

from enum import Enum, IntEnum


class MutationType(Enum):
    SYNONYMOUS = 'syn'
    MISSENSE = 'mis'
    NONSENSE = 'non'


class VariantType(IntEnum):
    INSERTION = 0
    DELETION = 1
    SUBSTITUTION = 2
    UNKNOWN = 3


class VariantClassification(IntEnum):
    CLASSIFIED = 0
    UNCLASSIFIED = 1
    MONOMORPHIC = 2


class ExecMode(str, Enum):
    SGE = 'sge'
    CDNA = 'cdna'


class SrcType(str, Enum):
    REF = 'ref'
    CDNA = 'cdna'
