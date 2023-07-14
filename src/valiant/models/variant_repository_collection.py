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
import logging
import sys
from typing import FrozenSet, Optional

from pyranges import PyRanges

from .pam_protection import PamProtectionVariantRepository
from .refseq_ranges import ReferenceSequenceRangeCollection
from .sge_config import SGEConfig
from .variant import VariantRepository


def _load_pam_protection_vcf(sgrna_ids: FrozenSet[str], pam: Optional[str]) -> PamProtectionVariantRepository:
    logging.debug("sgRNA ID's: %s" % ', '.join(sorted(sgrna_ids)) if sgrna_ids else "No sgRNA ID's.")

    vr = PamProtectionVariantRepository(sgrna_ids=sgrna_ids)

    if pam:
        if len(sgrna_ids) > 0:
            try:
                vr.load(pam)
            except ValueError as ex:
                logging.critical(ex.args[0])
                logging.critical("Failed to load the PAM protection variants!")
                sys.exit(1)
        else:
            logging.warning(
                "PAM protection variant file ignored: "
                "no sgRNA ID's associated to any targeton!")
    return vr


def _load_variant_repository(vcf_fp: Optional[str], ref_ranges: PyRanges, label: str, from_manifest: bool = False) -> Optional[VariantRepository]:
    if not vcf_fp:
        return None

    logging.debug(f"Loading {label} variants...")
    try:
        return (
            VariantRepository.load if from_manifest else
            VariantRepository.load_vcf
        )(vcf_fp, ref_ranges)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical(f"Failed to load {label} variants!")
        sys.exit(1)
    except FileNotFoundError as ex:
        logging.critical(ex.args[0])
        logging.critical(f"Failed to load {label} variants!")
        sys.exit(1)


@dataclass
class VariantRepositoryCollection:
    pam: PamProtectionVariantRepository
    background: Optional[VariantRepository]
    custom: Optional[VariantRepository]

    @classmethod
    def load(cls, config: SGEConfig, rsrs: ReferenceSequenceRangeCollection):

        # Load PAM protection variants
        pam_repository = _load_pam_protection_vcf(rsrs.sgrna_ids, config.pam_fp)

        # Load background variants
        background_variant_repository = (
            _load_variant_repository(
                config.bg_fp,
                rsrs._ref_ranges,
                'background'
            )
        )

        # Load custom variants
        variant_repository = (
            _load_variant_repository(
                config.vcf_fp,
                rsrs._ref_ranges,
                'custom',
                from_manifest=True
            )
        )

        return cls(pam_repository, background_variant_repository, variant_repository)
