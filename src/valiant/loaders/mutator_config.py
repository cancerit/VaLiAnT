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

from __future__ import annotations

from dataclasses import dataclass
import re

from .errors import InvalidMutator
from ..mutator_type import MutatorType
from ..int_pattern_builder import IntPatternBuilder
from ..utils import parse_opt_int_group


MUTATOR_TYPE_RE = {
    MutatorType.DEL: re.compile(r"^(\d+)del(\d*)$")
}


@dataclass(frozen=True, slots=True)
class MutatorConfig:
    type: MutatorType
    pt: IntPatternBuilder | None

    @classmethod
    def parse(cls, s: str) -> MutatorConfig:
        for t, pattern in MUTATOR_TYPE_RE.items():
            m = pattern.match(s)
            if m:

                # Parametric pattern
                try:
                    span = int(m.group(1))
                except ValueError:
                    raise InvalidMutator(f"Invalid mutator span in '{s}'!")

                offset = parse_opt_int_group(m, 2)

                try:
                    return cls(t, IntPatternBuilder(offset, span))
                except ValueError:
                    raise InvalidMutator(f"Invalid mutator parameters in '{s}'!")

        try:

            # Fixed pattern
            return cls(MutatorType(s), None)

        except ValueError:
            raise InvalidMutator(f"Invalid mutator type '{s}'!")
