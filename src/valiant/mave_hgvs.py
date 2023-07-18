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

"""
MAVE-HGVS mutation names

Limitations:

- liminal insertions are not recognised as exceptions (see HGVS 3' rule)

Ref.:
- https://www.mavedb.org/docs/mavehgvs/spec.html
- http://varnomen.hgvs.org
"""

from enum import Enum
from typing import NoReturn, Optional
from .enums import VariantType


class MAVEPrefix(Enum):
    LINEAR_GENOMIC = 'g'


# Encoded variant types
var_type_sub: int = VariantType.SUBSTITUTION.value
var_type_del: int = VariantType.DELETION.value
var_type_ins: int = VariantType.INSERTION.value


def _get_del_position(start: int, ref_length: int) -> str:
    return str(start) if ref_length == 1 else f"{start}_{start + ref_length - 1}"


def _get_delin_mave_nt_suffix(start: int, ref_length: int, alt: str) -> str:
    """Generate MAVE-HGVS name for a deletion-insertion"""

    return f"{_get_del_position(start, ref_length)}delins{alt}"


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
        _get_delin_mave_nt_suffix(start, ref_length, alt)
    )


# TODO: add special rules for 5' and 3' end insertions
def _get_insertion_mave_nt_suffix(start: int, alt: Optional[str]) -> str:
    """Generate MAVE-HGVS name for an insertion"""

    if alt is None or len(alt) == 0:
        raise ValueError("Invalid insertion: missing alternate!")

    return f"{start - 1}_{start}ins{alt}"


def _raise_invalid_deletion() -> NoReturn:
    raise ValueError("Invalid deletion: missing reference!")


def _get_deletion_mave_nt_suffix(start: int, ref: Optional[str], alt: Optional[str]) -> str:
    """Generate MAVE-HGVS name for a deletion"""

    if ref is None:
        _raise_invalid_deletion()

    ref_length: int = len(ref)
    n: int = len(ref)
    is_delin: bool = alt is not None

    if n == 0:
        _raise_invalid_deletion()

    return (
        _get_delin_mave_nt_suffix(start, ref_length, alt) if is_delin else
        _get_del_position(start, n) + 'del'
    )


def get_mave_nt(
    prefix: MAVEPrefix,
    var_type: int,
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
        _get_deletion_mave_nt_suffix(start, ref, alt) if var_type == var_type_del else
        None
    )

    if suffix is None:
        raise ValueError("Invalid variant type!")

    return f"{prefix.value}.{suffix}"
