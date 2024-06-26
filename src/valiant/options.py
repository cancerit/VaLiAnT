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

from dataclasses import dataclass

from .strings.strand import Strand


@dataclass(slots=True)
class Options:
    revcomp_minus_strand: bool
    include_no_op_oligo: bool
    oligo_max_length: int
    oligo_min_length: int
    allow_frame_shift: bool
    allow_non_syn: bool

    def should_rc(self, strand: Strand) -> bool:
        return self.revcomp_minus_strand and strand.is_minus
