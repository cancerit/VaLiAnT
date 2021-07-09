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

from enum import Enum


class TargetonMutator(Enum):
    DEL1 = '1del'
    DEL2_0 = '2del0'
    DEL2_1 = '2del1'
    IN_FRAME = 'inframe'
    SNV = 'snv'
    SNV_RE = 'snvre'
    STOP = 'stop'
    ALA = 'ala'
    AA = 'aa'


class MutationType(Enum):
    SYNONYMOUS = 0
    MISSENSE = 1
    NONSENSE = 2


class VariantType(Enum):
    INSERTION = 0
    DELETION = 1
    SUBSTITUTION = 2
    UNKNOWN = 3


class VariantClassification(Enum):
    CLASSIFIED = 0
    UNCLASSIFIED = 1
    MONOMORPHIC = 2
