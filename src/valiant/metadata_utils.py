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

import logging
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd

from .models.targeton import PamProtCDSTargeton
from .constants import META_MSEQ_NO_ADAPT_NO_RC, META_MUT_POSITION, META_NEW, META_PAM_CODON_ALT, META_PAM_CODON_MASK, META_PAM_CODON_REF, META_PAM_MUT_SGRNA_ID, META_PAM_MUT_START, META_PAM_SEQ, META_REF, META_REF_AA, META_REF_SEQ, META_REF_START, META_VCF_VAR_IN_CONST, METADATA_PAM_FIELDS
from .mave_hgvs import MAVEPrefix, get_mave_nt
from .utils import init_nullable_int_field, init_string_field


# Metadata table field used to determine whether an entry maps to a CDS
META_CDS_PROBE_FIELD = META_REF_AA


def is_metadata_row_cds(r: pd.Series) -> bool:
    return META_CDS_PROBE_FIELD in r and not pd.isna(r[META_CDS_PROBE_FIELD])


def _get_mave_nt_from_row(r: pd.Series, ref_start_field: str, ref_field: str, alt_field: str) -> str:
    """Generate MAVE-HGVS string with targeton-relative position"""

    return get_mave_nt(
        MAVEPrefix.LINEAR_GENOMIC,
        r.var_type,
        r[ref_start_field] - r.ref_start + 1,
        r[ref_field],
        r[alt_field])


def get_mave_nt_from_row(r: pd.Series) -> str:
    return _get_mave_nt_from_row(r, META_MUT_POSITION, META_REF, META_NEW)


def _get_mave_nt_pam_from_row(r: pd.Series) -> str:
    return _get_mave_nt_from_row(r, META_PAM_MUT_START, META_PAM_CODON_REF, META_PAM_CODON_ALT)


def get_mave_nt_pam_from_row(r: pd.Series) -> str:
    return (
        _get_mave_nt_pam_from_row if r[META_PAM_CODON_MASK] == 1 else
        get_mave_nt_from_row
    )(r)


def get_mave_nt_ref_from_row(r: pd.Series) -> str:
    return _get_mave_nt_from_row(r, META_MUT_POSITION, META_REF_NO_PAM, META_NEW)


def _get_field(
    df: pd.DataFrame,
    field: str,
    mask: Optional[pd.Series] = None
) -> pd.Series:
    return df.loc[mask, field] if mask is not None else df[field]


def get_string_lengths(df: pd.DataFrame, col_name: str, mask: Optional[pd.Series] = None) -> pd.Series:
    """Get string lengths"""

    return (
        _get_field(df, col_name, mask=mask)
        .astype('string')
        .str.len()
        .fillna(0)
        .astype(np.int32)
    )


def set_string_length_field(
    df: pd.DataFrame,
    string_col_name: str,
    length_col_name: str,
    mask: Optional[pd.Series] = None
) -> None:
    if mask is not None:
        df.loc[mask, length_col_name] = get_string_lengths(df, string_col_name, mask=mask)
    else:
        df[length_col_name] = get_string_lengths(df, string_col_name)


def get_slice(s: str, a: int, b: int) -> str:
    t = (
        s[a] if a == b else
        s[a:b + 1]
    )
    return t


def slice_string_field(
    df: pd.DataFrame,
    string_col_name: str,
    start_col_name: str,
    end_col_name: str,
    mask: Optional[pd.Series] = None
) -> List[str]:
    return [
        get_slice(s, a, b)
        for s, (a, b) in zip(
            _get_field(df, string_col_name, mask=mask),
            zip(
                _get_field(df, start_col_name, mask=mask),
                _get_field(df, end_col_name, mask=mask)
            ))
    ]


def set_slice_string_field(
    df: pd.DataFrame,
    string_col_name: str,
    start_col_name: str,
    end_col_name: str,
    slice_col_name: str,
    mask: Optional[pd.Series] = None
) -> None:
    if mask is not None:
        df.loc[mask, slice_col_name] = slice_string_field(
            df, string_col_name, start_col_name, end_col_name, mask=mask)
    else:
        df[slice_col_name] = slice_string_field(
            df, string_col_name, start_col_name, end_col_name)


# Temporary metadata fields
META_CDS_START = 'cds_start'
META_CDS_END = 'cds_end'
META_REF_END_POS = 'ref_end_pos'
META_REF_LENGTH = 'ref_length'
META_ALT_LENGTH = 'alt_length'
META_ALT_REF_DIFF = 'alt_ref_diff'
META_START_OFFSET = 'start_offset'
META_REF_END_OFFSET = 'ref_end_offset'
META_ALT_END_OFFSET = 'alt_end_offset'
META_REF_NO_PAM = 'ref_no_pam'

NO_CATEGORY: int = -1


def _init_pam_extended_required_fields(all_mutations: pd.DataFrame) -> None:
    rown: int = all_mutations.shape[0]

    # Initialise mask
    all_mutations[META_PAM_CODON_MASK] = np.zeros(rown, dtype=np.int8)

    # Initialise extended REF start positions
    init_nullable_int_field(all_mutations, META_PAM_MUT_START)

    # Initialise extended REF and ALT fields
    for col_name in [META_PAM_CODON_REF, META_PAM_CODON_ALT]:
        init_string_field(all_mutations, col_name)

    for field in METADATA_PAM_FIELDS:
        assert field in all_mutations


def _init_pam_extended_fields(all_mutations: pd.DataFrame) -> None:
    rown: int = all_mutations.shape[0]

    # Initialise start and end targetons, if any (-1 as default)
    all_mutations[META_CDS_START] = NO_CATEGORY * np.ones(rown, dtype=np.int8)
    all_mutations[META_CDS_END] = all_mutations[META_CDS_START].copy()

    # Initialise nullable integer fields
    for col_name in [
        META_ALT_REF_DIFF,
        META_ALT_END_OFFSET
    ]:
        init_nullable_int_field(all_mutations, col_name)

    # Initialise non-nullable integer field
    for col_name in [
        META_START_OFFSET,
        META_ALT_LENGTH
    ]:
        all_mutations[col_name] = np.zeros(rown, dtype=np.int32)


def get_interval_index(ts: List[Tuple[int, int]]) -> pd.IntervalIndex:
    return pd.IntervalIndex.from_tuples(ts, closed='both')


def _set_targeton_at_ref_start_end(
    all_mutations: pd.DataFrame,
    pam_codon_mask: pd.Series,
    pam_prot_cds_targetons: List[PamProtCDSTargeton]
) -> None:

    # Generate interval index to bin mutation positions
    pam_prot_cds_interval_index = pd.IntervalIndex.from_tuples([
        targeton.pos_range.to_tuple()
        for targeton in pam_prot_cds_targetons
    ], closed='both')

    # Assign REF start CDS (if any)
    logging.debug("Assigning targetons at REF starts...")
    all_mutations.loc[pam_codon_mask, META_CDS_START] = pd.cut(
        all_mutations.loc[pam_codon_mask, META_MUT_POSITION],
        pam_prot_cds_interval_index).cat.codes

    # Assign REF end CDS (if any)
    logging.debug("Assigning targetons at REF ends...")
    all_mutations.loc[pam_codon_mask, META_CDS_END] = pd.cut(
        all_mutations.loc[pam_codon_mask, META_REF_END_POS],
        pam_prot_cds_interval_index).cat.codes


def set_ref_meta(meta: pd.DataFrame) -> None:
    """
    Set temporary metadata fields describing REF

    Required fields:
    - META_REF: mutation reference
    - META_REF_START: reference sequence start genomic position
    - META_MUT_POSITION: mutation start genomic position

    Set fields:
    - META_REF_LENGTH: length of the mutation reference
    - META_REF_END_POS: mutation end genomic position
    - META_REF_END_OFFSET: mutation reference end relative position
    """

    # Initialise REF field (without PAM protection variants)
    init_string_field(meta, META_REF_NO_PAM)

    # Compute REF length
    set_string_length_field(meta, META_REF, META_REF_LENGTH)

    # Compute relative start position
    meta[META_START_OFFSET] = meta[META_MUT_POSITION].sub(meta[META_REF_START])

    # Compute REF end position
    # NOTE: reference lengths are corrected so that neither length zero or
    # length one REF's may incorrectly shift the inclusive end position.
    meta[META_REF_END_POS] = meta[META_MUT_POSITION].add(
        meta[META_REF_LENGTH].sub(1).clip(lower=0))

    meta[META_REF_END_OFFSET] = meta[META_REF_END_POS].sub(meta[META_REF_START])


def set_ref(meta: pd.DataFrame) -> None:
    """
    Set temporary metadata field with the REF before PAM protection

    Required fields:
    - META_START_OFFSET
    - META_REF_LENGTH
    - META_REF_SEQ
    - META_REF_END_OFFSET

    Set fields:
    - META_REF_NO_PAM
    """

    def set_slice(string_col_name: str, end_col_name: str, slice_col_name: str) -> None:
        set_slice_string_field(
            meta,
            string_col_name,
            META_START_OFFSET,
            end_col_name,
            slice_col_name,
            mask=meta[META_REF_LENGTH].gt(0))

    logging.debug("Assigning extended REF slice (no PAM protection)...")
    set_slice(
        META_REF_SEQ,
        META_REF_END_OFFSET,
        META_REF_NO_PAM)


def set_pam_extended_ref_alt(
    all_mutations: pd.DataFrame,
    pam_prot_cds_targetons: List[PamProtCDSTargeton]
) -> None:

    # Initialise required new fields
    _init_pam_extended_required_fields(all_mutations)


    # Filter for mutations that overlap at least one PAM-protected codon
    pam_codon_mask: pd.Series = all_mutations[META_PAM_MUT_SGRNA_ID].str.len() > 0

    if not pam_codon_mask.any():
        logging.warning("No sgRNA ID's assigned to PAM-protected CDS targetons!")
        return

    # Initialise extra fields
    _init_pam_extended_fields(all_mutations)

    def get_targeton(targeton_index: int) -> Optional[PamProtCDSTargeton]:
        assert targeton_index >= -1  # DEBUG
        return (
            pam_prot_cds_targetons[targeton_index] if targeton_index != NO_CATEGORY else
            None
        )

    def assign_pam_ref_alt_start(df: pd.DataFrame) -> pd.DataFrame:
        """Given a dataframe grouped by the start targeton, set the REF and ALT start"""

        targeton = get_targeton(df.name)

        if targeton:
            df[META_PAM_MUT_START] = df[META_MUT_POSITION].apply(
                lambda x: targeton.get_pam_ext_start(x)).astype(np.int32)

            df[META_START_OFFSET] = df[META_PAM_MUT_START].sub(df[META_REF_START])

        return df

    def assign_pam_ref_alt_end(df: pd.DataFrame) -> pd.Series:
        """Given a dataframe grouped by the end targeton, set the REF and ALT ends"""

        targeton = get_targeton(df.name)

        if targeton:
            df[META_REF_END_OFFSET] = df[META_REF_END_POS].apply(
                lambda x: targeton.get_pam_ext_ref_end(x)).sub(df[META_REF_START])

        df[META_ALT_END_OFFSET] = df[META_REF_END_OFFSET].add(df[META_ALT_REF_DIFF])

        return df

    def set_length_field(string_col_name: str, length_col_name: str) -> None:
        set_string_length_field(
            all_mutations, string_col_name, length_col_name, mask=pam_codon_mask)

    # Compute ALT length
    set_length_field(META_NEW, META_ALT_LENGTH)

    # Compute length difference between ALT and REF
    all_mutations.loc[pam_codon_mask, META_ALT_REF_DIFF] = (
        all_mutations.loc[pam_codon_mask, META_ALT_LENGTH] -
        all_mutations.loc[pam_codon_mask, META_REF_LENGTH]
    )

    _set_targeton_at_ref_start_end(
        all_mutations, pam_codon_mask, pam_prot_cds_targetons)

    # Assign extended REF and ALT start offsets
    all_mutations.loc[pam_codon_mask, :] = (
        all_mutations.loc[pam_codon_mask, :]
        .groupby(META_CDS_START, group_keys=False)
        .apply(assign_pam_ref_alt_start)
    )

    # Assign extended REF and ALT end offsets
    all_mutations.loc[pam_codon_mask, :] = (
        all_mutations.loc[pam_codon_mask, :]
        .groupby(META_CDS_END, group_keys=False)
        .apply(assign_pam_ref_alt_end)
    )

    def set_slice(string_col_name: str, end_col_name: str, slice_col_name: str) -> None:
        set_slice_string_field(
            all_mutations,
            string_col_name,
            META_START_OFFSET,
            end_col_name,
            slice_col_name,
            mask=pam_codon_mask)

    logging.debug("Assigning extended REF and ALT slices...")
    set_slice(
        META_PAM_SEQ,
        META_REF_END_OFFSET,
        META_PAM_CODON_REF)

    set_slice(
        META_MSEQ_NO_ADAPT_NO_RC,
        META_ALT_END_OFFSET,
        META_PAM_CODON_ALT)

    # Add mask (to filter when generating the VCF)
    all_mutations.loc[pam_codon_mask, META_PAM_CODON_MASK] = 1

    # Fill NA's in mask for custom variants in constant regions
    all_mutations[META_VCF_VAR_IN_CONST] = (
        all_mutations[META_VCF_VAR_IN_CONST].fillna(0) if META_VCF_VAR_IN_CONST in all_mutations else
        np.zeros(all_mutations.shape[0], dtype=pd.Int8Dtype)
    )
