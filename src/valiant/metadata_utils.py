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
from typing import List, Optional
import numpy as np
import pandas as pd

from .models.targeton import PamProtCDSTargeton
from .constants import META_MSEQ_NO_ADAPT_NO_RC, META_MUT_POSITION, META_MUTATOR, META_NEW, META_PAM_CODON_ALT, META_PAM_CODON_REF, META_PAM_MUT_SGRNA_ID, META_PAM_MUT_START, META_PAM_SEQ, META_REF, META_REF_AA, META_REF_START
from .mave_hgvs import MAVEPrefix, get_mave_nt


# Metadata table field used to determine whether an entry maps to a CDS
META_CDS_PROBE_FIELD = META_REF_AA


def is_metadata_row_cds(r: pd.Series) -> bool:
    return META_CDS_PROBE_FIELD in r and not pd.isna(r[META_CDS_PROBE_FIELD])


def get_mave_nt_from_row(r: pd.Series) -> str:
    return get_mave_nt(
        MAVEPrefix.LINEAR_GENOMIC,
        r.var_type,
        r.mut_position,
        r.ref,
        r.new)


def get_string_lengths(df: pd.DataFrame, mask: pd.Series, col_name: str) -> pd.Series:
    return (
        df.loc[mask, col_name]
        .astype('string')
        .str.len()
        .fillna(0)
        .astype(pd.Int32Dtype())
    )


def set_string_length_field(
    df: pd.DataFrame,
    mask: pd.Series,
    string_col_name: str,
    length_col_name: str
) -> None:
    df.loc[mask, length_col_name] = get_string_lengths(df, mask, string_col_name)


def get_slice(s: str, a: int, b: int) -> str:
    t = (
        s[a] if a == b else
        s[a:b + 1]
    )
    return t


def slice_string_field(
    df: pd.DataFrame,
    mask: pd.Series,
    string_col_name: str,
    start_col_name: str,
    end_col_name: str
) -> List[str]:
    return [
        get_slice(s, a, b)
        for s, (a, b) in zip(
            df.loc[mask, string_col_name],
            zip(
                df.loc[mask, start_col_name],
                df.loc[mask, end_col_name]
            ))
    ]


def set_slice_string_field(
    df: pd.DataFrame,
    mask: pd.Series,
    string_col_name: str,
    start_col_name: str,
    end_col_name: str,
    slice_col_name: str
) -> None:
    df.loc[mask, slice_col_name] = slice_string_field(
        df, mask, string_col_name, start_col_name, end_col_name)


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

NO_CATEGORY: int = -1


def _init_pam_extended_fields(all_mutations: pd.DataFrame) -> None:
    rown: int = all_mutations.shape[0]

    # Initialise start and end targetons, if any (-1 as default)
    all_mutations[META_CDS_START] = NO_CATEGORY * np.ones(rown, dtype=np.int8)
    all_mutations[META_CDS_END] = all_mutations[META_CDS_START].copy()

    # Initialise nullable integer fields
    for col_name in [
        META_PAM_MUT_START,
        META_REF_LENGTH,
        META_ALT_LENGTH,
        META_ALT_REF_DIFF,
        META_START_OFFSET,
        META_REF_END_OFFSET,
        META_ALT_END_OFFSET
    ]:
        all_mutations[col_name] = np.empty(rown, dtype=pd.Int32Dtype)

    # Initialise extended REF and ALT fields
    for col_name in [
        META_PAM_CODON_REF,
        META_PAM_CODON_ALT
    ]:
        all_mutations[col_name] = None
        all_mutations[col_name] = all_mutations[col_name].astype('string')


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


def set_pam_extended_ref_alt(
    all_mutations: pd.DataFrame,
    pam_prot_cds_targetons: List[PamProtCDSTargeton],
    debug_hash: Optional[str] = None
) -> None:

    # Initialise new fields
    _init_pam_extended_fields(all_mutations)

    # Filter for mutations that overlap at least one PAM-protected codon
    pam_codon_mask: pd.Series = all_mutations[META_PAM_MUT_SGRNA_ID].str.len() > 0

    def get_targeton(targeton_index: int) -> Optional[PamProtCDSTargeton]:
        assert targeton_index >= -1  # DEBUG
        return (
            pam_prot_cds_targetons[targeton_index] if targeton_index != NO_CATEGORY else
            None
        )

    def assign_pam_ref_alt_start(df: pd.DataFrame) -> pd.DataFrame:
        """Given a dataframe grouped by the start targeton, set the REF and ALT start"""

        targeton = get_targeton(df.name)

        df[META_PAM_MUT_START] = (
            df[META_MUT_POSITION].apply(
                lambda x: targeton.get_pam_ext_start(x)) if targeton else
            df[META_MUT_POSITION]
        )
        df[META_START_OFFSET] = df[META_PAM_MUT_START].sub(df[META_REF_START])

        return df

    def assign_pam_ref_alt_end(df: pd.DataFrame) -> pd.Series:
        """Given a dataframe grouped by the end targeton, set the REF and ALT ends"""

        targeton = get_targeton(df.name)

        df[META_REF_END_OFFSET] = (
            df[META_REF_END_POS].apply(
                lambda x: targeton.get_pam_ext_ref_end(x)) if targeton else
            df[META_REF_END_POS]
        ).sub(df[META_REF_START])

        df[META_ALT_END_OFFSET] = df[META_REF_END_OFFSET].add(df[META_ALT_REF_DIFF])

        return df

    def set_length_field(string_col_name: str, length_col_name: str) -> None:
        set_string_length_field(
            all_mutations, pam_codon_mask, string_col_name, length_col_name)

    # Compute REF length
    set_length_field(META_REF, META_REF_LENGTH)

    # Compute ALT length
    set_length_field(META_NEW, META_ALT_LENGTH)

    # Compute length difference between ALT and REF
    all_mutations.loc[pam_codon_mask, META_ALT_REF_DIFF] = (
        all_mutations.loc[pam_codon_mask, META_ALT_LENGTH] -
        all_mutations.loc[pam_codon_mask, META_REF_LENGTH]
    )

    # Compute REF end position
    # NOTE: refernce lengths are corrected so that neither length zero or
    # length one REF's may incorrectly shift the inclusive end position.
    all_mutations.loc[pam_codon_mask, META_REF_END_POS] = (
        all_mutations.loc[pam_codon_mask, META_MUT_POSITION] +
        all_mutations.loc[pam_codon_mask, META_REF_LENGTH].sub(1).clip(lower=0)
    )

    _set_targeton_at_ref_start_end(
        all_mutations, pam_codon_mask, pam_prot_cds_targetons)

    # Assign extended REF and ALT start offsets
    all_mutations.loc[pam_codon_mask,:] = (
        all_mutations.loc[pam_codon_mask,:]
        .groupby(META_CDS_START, group_keys=False)
        .apply(assign_pam_ref_alt_start)
    )

    # Assign extended REF and ALT end offsets
    all_mutations.loc[pam_codon_mask,:] = (
        all_mutations.loc[pam_codon_mask,:]
        .groupby(META_CDS_END, group_keys=False)
        .apply(assign_pam_ref_alt_end)
    )

    logging.debug("Assigning extended REF and ALT slices...")
    def set_slice(string_col_name: str, end_col_name: str, slice_col_name: str) -> None:
        set_slice_string_field(
            all_mutations,
            pam_codon_mask,
            string_col_name,
            META_START_OFFSET,
            end_col_name,
            slice_col_name)

    set_slice(
        META_PAM_SEQ,
        META_REF_END_OFFSET,
        META_PAM_CODON_REF)

    set_slice(
        META_MSEQ_NO_ADAPT_NO_RC,
        META_ALT_END_OFFSET,
        META_PAM_CODON_ALT)
