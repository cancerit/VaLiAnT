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

from typing import Dict, List, Optional, Set
import pandas as pd
from .constants import META_MUT_POSITION, META_REF, META_PAM_MUT_SGRNA_ID, ARRAY_SEPARATOR


def get_sgrna_ids(frame: int, codon_to_sgrna_ids: Dict[int, str], mut_pos: int, mut_ref: Optional[str]) -> Set[str]:
    """List the identifiers of the sgRNA's whose variants map to codons altered by the mutation"""

    start: int = mut_pos + frame
    codon_start: int = start // 3
    ref_length: int = len(mut_ref) if mut_ref else 0

    return (
        {
            codon_to_sgrna_ids[codon_index]
            for codon_index in range(codon_start, ((start + ref_length) // 3) + 1)
            if codon_index in codon_to_sgrna_ids
        } if ref_length > 1 else
        {codon_to_sgrna_ids[codon_start]} if codon_start in codon_to_sgrna_ids else
        set()
    )


def get_sgrna_ids_from_row(frame: int, codon_to_sgrna_ids: Dict[int, str], r: pd.Series) -> Set[str]:
    """
    List the identifiers of the sgRNA's whose variants map to codons altered by the mutation

    Metadata table row expected, with relative start position and reference sequence.
    """

    mut_ref: Optional[str] = r[META_REF] if not pd.isna(r[META_REF]) else None
    return get_sgrna_ids(frame, codon_to_sgrna_ids, r[META_MUT_POSITION], mut_ref)


def sgrna_ids_to_string(sgrna_ids: Set[str]) -> str:
    return ARRAY_SEPARATOR.join(sorted(sgrna_ids))


def get_sgrna_ids_from_row_as_string(frame: int, codon_to_sgrna_ids: Dict[int, str], r: pd.Series) -> str:
    """Concatenated sgRNA ID's or empty string"""

    return sgrna_ids_to_string(
        get_sgrna_ids_from_row(
            frame, codon_to_sgrna_ids, r))


def set_metadata_sgrna_ids(frame: int, codon_to_sgrna_ids: Dict[int, str], metadata: pd.DataFrame) -> None:
    """Set the sgRNA ID field in the metadata table (SGE only)"""

    def get_sgrna_id_str(r: pd.Series) -> str:
        return get_sgrna_ids_from_row_as_string(
            frame, codon_to_sgrna_ids, r)

    metadata[META_PAM_MUT_SGRNA_ID] = metadata.apply(
        get_sgrna_id_str, axis=1).astype('string')


def get_codon_index(seq_frame: int, seq_start: int, x: int) -> int:
    return (x + seq_frame - seq_start) // 3


def get_codon_indices_in_range(
    seq_frame: int,
    seq_start: int,
    seq_end: int,
    target_start: int,
    target_end: int
) -> List[int]:
    """Get the indices of the codons spanned by the target range, if any"""

    def _get_codon_index(x: int) -> int:
        return get_codon_index(seq_frame, seq_start, x)

    # NOTE: in a `PositionRange`, `end` is guaranteed to be larger than `start`
    assert seq_end >= seq_start and target_end >= target_start

    if target_end < seq_start or target_start > seq_end:
        return []

    # Get first codon index
    codon_start: int = 0 if target_start <= seq_start else _get_codon_index(target_start)

    if target_end == target_start:
        return [codon_start]

    # Get last codon index
    codon_end: int = (
        _get_codon_index(seq_end) if target_end >= seq_end else
        _get_codon_index(target_end)
    )

    # Generate full codon index range
    return (
        [codon_start] if codon_end == codon_start else
        list(range(codon_start, codon_end + 1))
    )
