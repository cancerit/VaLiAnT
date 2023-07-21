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
from dataclasses import dataclass
import logging
from typing import Dict, List, Optional, Set, FrozenSet
import pandas as pd
from pyranges import PyRanges

from ..loaders.vcf import get_vcf
from ..utils import get_id_column
from .variant import get_variant, SubstitutionVariant


@dataclass(frozen=True)
class PamVariant(SubstitutionVariant):
    __slots__ = {'genomic_position', 'ref', 'alt', 'sgrna_id'}

    sgrna_id: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if len(self.ref) != 1 or len(self.alt) != 1:
            raise ValueError(
                "Only single-nucleotide substitutions are allowed "
                "for the purposes of PAM protection!")

    @classmethod
    def from_substitution(cls, var: SubstitutionVariant, sgrna_id: str) -> PamVariant:
        if not isinstance(var, SubstitutionVariant):
            raise ValueError("PAM protection variants must be substitutions!")
        return PamVariant(var.genomic_position, var.ref, var.alt, sgrna_id)


class PamProtectionVariantRepository:
    __slots__ = {'sgrna_ids', '_variants', '_ranges'}

    def __init__(self, sgrna_ids: FrozenSet[str] = None) -> None:
        self.sgrna_ids = sgrna_ids or frozenset()
        self._variants: Dict[str, Set[PamVariant]] = {
            sgrna_id: set()
            for sgrna_id in sgrna_ids
        } if sgrna_ids else {}
        self._ranges: Optional[PyRanges] = None

    @property
    def count(self) -> int:
        return sum(map(len, self._variants.values()))

    def get_sgrna_variants(self, sgrna_id: str) -> FrozenSet[PamVariant]:
        try:
            return frozenset(self._variants[sgrna_id])

        except KeyError:
            raise RuntimeError(f"sgRNA ID '{sgrna_id} not found!'")

    def get_sgrna_variants_bulk(self, sgrna_ids: FrozenSet[str]) -> FrozenSet[PamVariant]:
        empty: FrozenSet[PamVariant] = frozenset()
        return empty.union(*[
            self.get_sgrna_variants(sgrna_id)
            for sgrna_id in sgrna_ids
        ])

    def load(self, fp: str) -> None:

        # Load variants from VCF
        with get_vcf(fp) as variant_file:
            for record in variant_file.fetch():
                sgrna_id: str = record.info['SGRNA'].strip()
                if sgrna_id in self._variants:
                    self._variants[sgrna_id].add(
                        PamVariant.from_substitution(get_variant(record), sgrna_id))

        # Log loaded variant statistics
        logging.debug("Collected %d PAM protection variants." % self.count)
        for sgrna_id in self.sgrna_ids:
            if not self._variants[sgrna_id]:
                logging.info("No PAM protection variants for sgRNA %s." % sgrna_id)

        # Populate genomic range table
        df: pd.DataFrame = pd.DataFrame.from_records([
            (*variant.get_pyrange_record(), sgrna_id)
            for sgrna_id, variants in self._variants.items()
            for variant in variants
        ], columns=[
            'Chromosome',
            'Start',
            'End',
            'sgrna_id'
        ])
        df.sgrna_id = df.sgrna_id.astype('category')
        df['variant_id'] = get_id_column(df.shape[0])
        self._ranges = PyRanges(df)


def get_position_to_sgrna_ids(pam_variants: List[PamVariant]) -> Dict[int, str]:
    return {
        variant.genomic_position.position: variant.sgrna_id
        for variant in pam_variants
    }
