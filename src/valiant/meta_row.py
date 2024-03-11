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

from dataclasses import dataclass

from .strings.dna_str import DnaStr
from .uint_range import UIntRange
from .variant import Variant


sql_select_meta = """
select
    original_start,
    ref_start,
    ref_end,
    ref,
    alt,
    ref_aa,
    alt_aa,
    vcf_var_id,
    vcf_alias,
    vcf_nt,
    mutator,
    in_const,
    oligo,
    mutation_type,
    start_exon_index,
    start_codon_index,
    start_ppe_start,
    end_exon_index,
    end_codon_index,
    end_ppe_start,
    sgrna_ids
from v_meta
order by original_start, ref_start, ref_end, ref, alt, vcf_var_id, mutator
"""


@dataclass(slots=True)
class MetaRow:
    ref_pos: int
    alt_pos: int
    end: int
    ref: str
    alt: str
    ref_aa: str
    alt_aa: str
    vcf_var_id: str
    vcf_alias: str
    vcf_nt: str
    mutator: str
    in_const: int
    oligo: str
    mutation_type: str
    start_exon_index: int | None
    start_codon_index: int | None
    start_ppe_start: int | None
    end_exon_index: int | None
    end_codon_index: int | None
    end_ppe_start: int | None
    sgrna_ids: str

    def _to_variant(self, pos: int) -> Variant:
        return Variant(pos, DnaStr(self.ref), DnaStr(self.alt))

    def to_ref_variant(self) -> Variant:
        return self._to_variant(self.ref_pos)

    def to_alt_variant(self) -> Variant:
        return self._to_variant(self.alt_pos)

    @property
    def alt_ref_range(self) -> UIntRange:
        return UIntRange(self.alt_pos, self.end)

    @property
    def overlaps_codon(self) -> bool:
        return (
            self.start_exon_index is not None or
            self.end_exon_index is not None
        )

    @property
    def pam_ref_start(self) -> int:
        if self.start_ppe_start is None:
            return self.alt_pos
        return min(self.alt_pos, self.start_ppe_start)

    @property
    def pam_ref_end(self) -> int:
        if self.end_ppe_start is None:
            return self.end
        return max(self.end, self.end_ppe_start)

    @property
    def pam_ref_range(self) -> UIntRange:
        return UIntRange(self.pam_ref_start, self.pam_ref_end)
