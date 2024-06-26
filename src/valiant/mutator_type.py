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

from enum import Enum


class MutatorType(str, Enum):
    # Parametric span
    DEL = 'del'
    # Single-nucleotide span
    SNV = 'snv'
    # Codon span
    SNV_RE = 'snvre'
    IN_FRAME = 'inframe'
    ALA = 'ala'
    STOP = 'stop'
    AA = 'aa'


PARAMETRIC_MUTATOR_TYPES: set[MutatorType] = {
    MutatorType.DEL
}


CDS_MUTATOR_TYPES: set[MutatorType] = {
    MutatorType.SNV_RE,
    MutatorType.ALA,
    MutatorType.STOP,
    MutatorType.AA,
    MutatorType.IN_FRAME
}


ANNOTABLE_MUTATOR_TYPES: set[MutatorType] = {
    MutatorType.SNV,
    MutatorType.SNV_RE,
    MutatorType.ALA,
    MutatorType.STOP,
    MutatorType.AA
}


# Mutator types based on the output of other mutators
DEPENDENT_MUTATOR_TYPES: dict[MutatorType, set[MutatorType]] = {
    MutatorType.SNV_RE: {MutatorType.SNV}
}
