########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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

"""
MAVE-HGVS mutation names

Ref.: https://www.mavedb.org/docs/mavehgvs/spec.html#sequence-prefixes-and-sequence-types
"""

from enum import Enum
from typing import Optional
from .enums import VariantType


class MAVEPrefix(Enum):
    LINEAR_GENOMIC = 'g'


# Encoded variant types
var_type_sub: int = VariantType.SUBSTITUTION.value
var_type_del: int = VariantType.DELETION.value
var_type_ins: int = VariantType.INSERTION.value


def _get_delin_mave_nt_suffix(start: int, alt: str) -> str:
    """Generate MAVE-HGVS name for a deletion-insertion"""

    return f"{start}delins{alt}"


def _get_snv_mave_nt_suffix(start: int, ref: str, alt: str) -> str:
    """Generate MAVE-HGVS name for a deletion-insertion"""

    return f"{start}{ref}>{alt}"


def _get_substitution_mave_nt_suffix(start: int, ref: Optional[str], alt: Optional[str]) -> str:
    ref_length: int = len(ref) if ref is not None else 0
    alt_length: int = len(alt) if alt is not None else 0

    if ref_length == 0:
        raise ValueError("Invalid substitution: missing reference!")
    if alt_length == 0:
        raise ValueError("Invalid substitution: missing alternate!")

    is_snv: bool = ref_length == 1 and alt_length == 1

    return (
        _get_snv_mave_nt_suffix(start, ref, alt) if is_snv else
        _get_delin_mave_nt_suffix(start, alt)
    )


def _get_insertion_mave_nt_suffix(start: int, alt: Optional[str]) -> str:
    """Generate MAVE-HGVS name for an insertion"""

    if alt is None or len(alt) == 0:
        raise ValueError("Invalid insertion: missing alternate!")
    return f"{start}ins{alt}"


def _get_deletion_mave_nt_suffix(start: int) -> str:
    """Generate MAVE-HGVS name for a deletion"""

    # TODO: verify format for multi-nucleotide deletions...
    return f"{start}del"


def get_mave_nt(
    prefix: MAVEPrefix,
    var_type: int,
    source: str,
    start: int,
    ref: Optional[str],
    alt: Optional[str]
) -> str:
    """
    Generate MAVE-HGVS name
    """

    if start < 0:
        raise ValueError(f"Invalid start position: {start}!")

    suffix: Optional[str] = (
        _get_substitution_mave_nt_suffix(start, ref, alt) if var_type == var_type_sub else
        _get_insertion_mave_nt_suffix(start, alt) if var_type == var_type_ins else
        _get_deletion_mave_nt_suffix(start) if var_type == var_type_del else
        None
    )

    if suffix is None:
        raise ValueError("Invalid variant type!")

    return f"{prefix.value}.{suffix}"
