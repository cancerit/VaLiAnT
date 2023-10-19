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

from __future__ import annotations

import logging
from dataclasses import dataclass

from pysam import VariantRecord

from .contig_filter import ContigFilter
from .enums import VariantType, VariantClassification
from .loaders.vcf import load_vcf
from .strings.dna_str import DnaStr
from .strings.nucleotide import Nucleotide
from .variant import Variant
from .vcf_variant import VcfVariant


def validate_vcf_record_ref(var_id: str | None, r: VariantRecord) -> None:
    if r.ref == '.':
        id_: str = f" {var_id}" if var_id is not None else ''
        raise ValueError(f"Empty REF for variant{id_} at {r.contig}:{r.pos}!")


@dataclass
class CustomVariant(Variant):
    contig: str
    id: str | None
    vcf_nt: str | None
    var_type: VariantType
    var_class: VariantClassification

    def as_variant(self) -> Variant:
        return Variant(self.pos, self.ref, self.alt)

    @classmethod
    def from_variant(
        cls,
        v: Variant,
        contig: str,
        id: str | None,
        vcf_nt: str | None,
        var_type: VariantType,
        var_class: VariantClassification
    ) -> CustomVariant:
        return cls(
            v.pos, v.ref, v.alt,
            contig, id, vcf_nt,
            var_type, var_class)

    @classmethod
    def from_record_with_id(cls, id: str | None, r: VariantRecord) -> CustomVariant:
        validate_vcf_record_ref(id, r)

        contig: str = r.contig
        start: int = r.pos
        ref = DnaStr(r.ref.upper())

        # Monomorphic variant
        if r.alts is None:
            return cls(
                start, ref, DnaStr.empty(), contig, id, None,
                VariantType.UNKNOWN, VariantClassification.MONOMORPHIC)

        # TODO: handle multiple ALT's?
        raw_var = VcfVariant(start, ref, DnaStr(r.alts[0].upper()))

        def classify(
            var: Variant,
            var_type: VariantType,
            var_class: VariantClassification = VariantClassification.CLASSIFIED,
            vcf_nt: Nucleotide | None = None
        ) -> CustomVariant:
            return cls.from_variant(var, contig, id, vcf_nt, var_type, var_class)

        alt_ref_delta: int = raw_var.alt_ref_delta

        if alt_ref_delta == 0:
            return classify(raw_var, VariantType.SUBSTITUTION)

        norm_var, nt = raw_var.normalise()

        if nt:
            return classify(norm_var, (
                VariantType.DELETION if alt_ref_delta < 0 else
                VariantType.INSERTION
            ), vcf_nt=nt)

        return classify(
            norm_var,
            VariantType.SUBSTITUTION,
            var_class=VariantClassification.UNCLASSIFIED)

    def _vcf_filter(self) -> bool:
        if self.var_class == VariantClassification.MONOMORPHIC:
            logging.info("Monomorphic variant at %s:%d (SKIPPED)." % (self.contig, self.pos))
            return False
        return True

    @classmethod
    def load_vcf(cls, fp: str, ft: ContigFilter, vcf_id_tag: str | None = None) -> list[CustomVariant]:
        return load_vcf(
            fp, ft, cls.from_record_with_id,
            filter_f=cls._vcf_filter,
            vcf_id_tag=vcf_id_tag)
