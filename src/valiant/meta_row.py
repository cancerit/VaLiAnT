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

from .enums import VariantType


sql_select_meta = """
select
    ref_start,
    ref,
    alt,
    ref_aa,
    alt_aa,
    vcf_var_id,
    vcf_alias,
    mutator,
    in_const,
    oligo,
    mutation_type,
    start_exon_index,
    start_codon_index,
    end_exon_index,
    end_codon_index,
    sgrna_ids
from v_meta
"""


@dataclass(slots=True)
class MetaRow:
    pos: int
    ref: str
    alt: str
    ref_aa: str
    alt_aa: str
    vcf_var_id: str
    vcf_alias: str
    mutator: str
    in_const: int
    oligo: str
    mutation_type: str
    start_exon_index: int | None
    start_codon_index: int | None
    end_exon_index: int | None
    end_codon_index: int | None
    sgrna_ids: str

    @property
    def variant_type(self) -> VariantType:
        has_ref = bool(self.ref)
        has_alt = bool(self.alt)
        if has_ref:
            return (
                VariantType.SUBSTITUTION if has_alt else
                VariantType.DELETION
            )
        if has_alt:
            return VariantType.INSERTION
        raise ValueError("Invalid metadata row!")
