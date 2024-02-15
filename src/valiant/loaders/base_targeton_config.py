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

import abc
from dataclasses import dataclass
from typing import TypeVar

from ..uint_range import UIntRange
from .mutator_config import MutatorConfig
from .utils import parse_list


def parse_mutators(s: str) -> list[MutatorConfig]:
    mutator_codes = sorted(set(parse_list(s)))
    return list(map(MutatorConfig.parse, mutator_codes))


@dataclass
class BaseTargetonConfig(abc.ABC):
    ref: UIntRange
    region_2: UIntRange

    @property
    @abc.abstractmethod
    def name(self) -> str:
        pass


BaseTargetonConfigT = TypeVar('BaseTargetonConfigT', bound=BaseTargetonConfig)
