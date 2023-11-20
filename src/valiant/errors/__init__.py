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

from abc import ABC, abstractmethod


class CustomException(Exception, ABC):

    @property
    @abstractmethod
    def message(self) -> str:
        pass


class NonUniqueError(Exception):
    pass


class SequenceNotFound(Exception):
    pass


class InvalidMutatorForTarget(ValueError):
    pass


class InvalidConfig(Exception):
    pass


class InvalidBackground(Exception):
    pass


class InvalidVariantRef(Exception):
    pass


class GenomicRangeOutOfBounds(Exception):
    pass