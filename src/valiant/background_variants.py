########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from .variant import RegisteredVariant


class InvalidBackgroundVariant(Exception):
    pass


@dataclass(slots=True)
class RegisteredBackgroundVariant(RegisteredVariant):
    start_exon_index: int | None
    start_codon_index: int | None
    end_exon_index: int | None
    end_codon_index: int | None

    @property
    def in_cds(self) -> bool:
        return (
            self.start_exon_index is not None or
            self.end_exon_index is not None
        )

    def get_non_syn_warn_message(self, contig, is_frame_shifting: bool) -> str:
        prefix = "Frame shifting" if is_frame_shifting else "Non-synonymous"
        return f"{prefix} background variant at {contig}:{self.pos}!"
