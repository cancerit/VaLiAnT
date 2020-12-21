########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
# 
# The Software remains the property of Genome Research Ltd ("GRL").
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
# 
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
# 
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
#############################

from typing import Callable, Dict, List, Tuple
import numpy as np
import pandas as pd
from .constants import DNA_ALPHABET_SET
from .enums import MutationType, VariantType
from .globals import NUCLEOTIDES, TRIPLET_RCS, TRIPLETS
from .models.base import GenomicRange
from .models.codon_table import CodonTable, STOP_CODE
from .string_mutators import replace_single_nucleotides
from .utils import get_inner_cds_relative_boundaries, validate_strand

SYN_CODON_FIELDS = [
    'triplet',
    'syn_triplet',
    'ref_aa',
    'codon_rank'
]

SNV_SYN_FIELDS = [
    'pos',
    'alt_triplet',
    'ref_aa',
    'alt_aa',
    'mut_type'
]

syn_code: int = MutationType.SYNONYMOUS.value
mis_code: int = MutationType.MISSENSE.value
non_code: int = MutationType.NONSENSE.value


def as_nucleotides(values) -> pd.Categorical:
    return pd.Categorical(values, categories=NUCLEOTIDES)


def as_triplets(values) -> pd.Categorical:
    return pd.Categorical(values, categories=TRIPLETS)


def as_amino_acids(values, amino_acid_symbols: List[str]) -> pd.Categorical:
    return pd.Categorical(values, categories=amino_acid_symbols)


def _get_triplet_aa_table(
    pairs: List[Tuple[str, str]],
    amino_acid_symbols: List[str],
    suffix: str
) -> pd.DataFrame:
    triplet_col: str = f"triplet_{suffix}"
    aa_col: str = f"aa_{suffix}"
    df: pd.DataFrame = pd.DataFrame.from_records(
        pairs, columns=[triplet_col, aa_col])
    df[triplet_col] = as_triplets(df[triplet_col])
    df[aa_col] = as_amino_acids(df[aa_col], amino_acid_symbols)
    return df


def build_aa_sub_table(codon_table: CodonTable, strand: str) -> pd.DataFrame:
    plus_strand: bool = strand == '+'
    amino_acid_symbols: List[str] = codon_table.amino_acid_symbols
    aa_count: int = len(amino_acid_symbols) - 1
    tr_table_f: Callable[[], List[Tuple[str, str]]] = (
        codon_table.get_translation_table if plus_strand else
        codon_table.get_translation_table_rc
    )
    tr_top_table_f: Callable[[], List[Tuple[str, str]]] = (
        codon_table.get_top_ranking_codon_table if plus_strand else
        codon_table.get_top_ranking_codon_table_rc
    )

    df: pd.DataFrame = _get_triplet_aa_table([
        (t, aa)
        for t, aa in tr_table_f()
        for _ in range(aa_count)
    ], amino_acid_symbols, 'ref')

    aa_df: pd.DataFrame = _get_triplet_aa_table([
        (t, aa)
        for aa, t in tr_top_table_f()
        if aa != STOP_CODE
    ], amino_acid_symbols, 'alt')

    df = pd.concat([
        df,
        pd.concat([aa_df] * len(TRIPLETS), axis=0, ignore_index=True)
    ], axis=1)

    return df[df.aa_ref != df.aa_alt].set_index('triplet_ref').rename(columns={
        'triplet_alt': 'alt',
        'aa_ref': 'ref_aa',
        'aa_alt': 'alt_aa'
    })


def get_fast_aa_subs(aa_sub_table: pd.DataFrame, frame: int, seq: str) -> pd.DataFrame:
    seq_len: int = len(seq)

    # Get longest in-frame CDS subsequence
    start, end = get_inner_cds_relative_boundaries(seq_len, frame)

    # Obtain reference triplets and relative positions
    df: pd.DataFrame = pd.DataFrame.from_records(
        [(seq[i:i + 3], i) for i in range(start, 3 * (end // 3), 3)],
        columns=['triplet', 'pos'])
    df.triplet = as_triplets(df.triplet)
    df.pos = df.pos.astype(np.int32)

    # Join alternative triplets and amino acid change information
    df = df.set_index('triplet', drop=True).join(aa_sub_table).reset_index()

    # Render mutated sequences
    df['mseq'] = _get_snvre_mseq(seq, df)

    # Add variant type
    df['var_type'] = np.int8(VariantType.SUBSTITUTION.value)

    return df.rename(columns={
        'triplet': 'ref',
        'pos': 'mut_position',
        'alt': 'new'
    })


def _build_base_snv_table() -> pd.DataFrame:
    data = pd.DataFrame.from_records([
        (t, *res)
        for t in TRIPLETS
        for res in replace_single_nucleotides(t, DNA_ALPHABET_SET)
    ], columns=['triplet', 'pos', 'ref', 'alt', 'alt_triplet'])

    # Compress data
    data.pos = data.pos.astype(np.int32)
    data.triplet = as_triplets(data.triplet)
    data.alt_triplet = as_triplets(data.alt_triplet)
    data.ref = as_nucleotides(data.ref)
    data.alt = as_nucleotides(data.alt)

    # Set index
    return data.set_index('triplet', drop=False)


def build_snv_table(base_snv_table: pd.DataFrame, codon_table: CodonTable, strand: str) -> pd.DataFrame:
    validate_strand(strand)

    tr: Callable[[str], str] = (
        codon_table.translate if strand == '+' else
        codon_table.translate_rc
    )

    df: pd.DataFrame = base_snv_table.copy()

    # Get amino acid information
    amino_acid_symbols: List[str] = codon_table.amino_acid_symbols
    df['ref_aa'] = as_amino_acids(df.triplet.apply(tr), amino_acid_symbols)
    df['alt_aa'] = as_amino_acids(df.alt_triplet.apply(tr), amino_acid_symbols)

    # Classify mutations
    df['mut_type'] = np.select([
        df.alt_aa == df.ref_aa,
        df.alt_aa == STOP_CODE
    ], [
        MutationType.SYNONYMOUS.value,
        MutationType.NONSENSE.value
    ], default=MutationType.MISSENSE.value).astype(np.int8)

    return df


def get_fast_exon_snv(
    snv_table: pd.DataFrame,
    cds_seq: str,
    genomic_range: GenomicRange,
    cds_prefix_length: int,
    cds_suffix_length: int,
    reset_index: bool = True
) -> pd.DataFrame:
    """
    Get pre-computed metadata for an in-frame coding sequence

    cds_seq: CDS-extended target sequence
    genomic_range: target start position
    """

    cds_seq_len: int = len(cds_seq)

    if cds_seq_len == 0:
        return snv_table.loc[[]]

    if cds_seq_len % 3 != 0:
        raise ValueError("The sequence length must be a multiple of three!")

    if cds_seq_len != len(genomic_range) + cds_prefix_length + cds_suffix_length:
        raise ValueError("Inconsistent CDS sequence and boundaries!")

    # Assemble SNV table
    df: pd.DataFrame = snv_table.loc[[
        cds_seq[i:i + 3]
        for i in range(0, cds_seq_len, 3)
    ]]

    # Adjust SNV positions
    df.pos += (
        np.arange(0, df.shape[0] // 3, step=3, dtype=np.int32).repeat(9)
        + np.int32(genomic_range.start - cds_prefix_length)
    )

    # Filter out SNV's exceeding the target range
    if cds_prefix_length > 0:
        df = df.loc[df.pos >= np.int32(genomic_range.start)]
    if cds_suffix_length > 0:
        df = df.loc[df.pos <= np.int32(genomic_range.end)]

    return df.reset_index(drop=True) if reset_index else df


def _get_base_all_syn_codons_plus(codon_table: CodonTable) -> pd.DataFrame:
    return pd.DataFrame.from_records([
        (triplet, syn_triplet, codon_table._codon2aa[triplet], rank)
        for triplet, syn_triplets in codon_table._codon2syn.items()
        for rank, syn_triplet in enumerate(syn_triplets)
    ], columns=SYN_CODON_FIELDS)


def _get_base_all_syn_codons_minus(codon_table: CodonTable) -> pd.DataFrame:
    """
    triplet (+), syn_rc_triplet_rc (+), syn_rc_triplet_aa (-), rank
    """
    return pd.DataFrame.from_records([
        (triplet, TRIPLET_RCS[syn_triplet], codon_table._codon2aa[syn_triplet], rank)
        for triplet, syn_triplets in codon_table._codonrc2syn.items()
        for rank, syn_triplet in enumerate(syn_triplets)
    ], columns=SYN_CODON_FIELDS)


STRAND_ALL_SYN_CODONS: Dict[str, Callable[[CodonTable], pd.DataFrame]] = {
    '+': _get_base_all_syn_codons_plus,
    '-': _get_base_all_syn_codons_minus
}


def get_all_syn_codons(codon_table: CodonTable, strand: str) -> pd.DataFrame:
    validate_strand(strand)

    df = STRAND_ALL_SYN_CODONS[strand](codon_table)

    # Compress data
    df.triplet = as_triplets(df.triplet)
    df.syn_triplet = as_triplets(df.syn_triplet)
    df.ref_aa = as_amino_acids(df.ref_aa, codon_table.amino_acid_symbols)
    df['alt_aa'] = df.ref_aa
    df.codon_rank = df.codon_rank.astype(np.int8)

    # Index by triplet
    return df.set_index('triplet')


def get_top_syn_codons(syn_df: pd.DataFrame) -> pd.DataFrame:
    return syn_df[syn_df.codon_rank == 0].drop('codon_rank', axis=1)


def _get_syn_from_snvs(syn_df: pd.DataFrame, snvs: pd.DataFrame, mut_type_code: int) -> pd.DataFrame:
    """
    Given a pre-computed metadata table for all or a subset of synonymous codons,
    join the SNV metadata (only with the fields listed in SNV_SYN_FIELDS) for a given
    mutation type with it.
    """

    # Filter SNV metadata by mutation type
    snvs_mut_type: pd.DataFrame = snvs[snvs.mut_type == mut_type_code]

    if snvs_mut_type.shape[0] == 0:
        return pd.DataFrame()

    # Join SNV mutations with the synonymous table
    df: pd.DataFrame = snvs_mut_type.join(
        syn_df.syn_triplet, on=['alt_triplet'])

    # Drop synonyms that are equal to the reference
    df = df[df.syn_triplet != df.index]

    # Drop variants with no possible synonymous variants
    df = df[~pd.isna(df.syn_triplet)]

    df['ref'] = df.index.values
    df['alt'] = df.syn_triplet
    return df.drop(['alt_triplet', 'syn_triplet'], axis=1).drop_duplicates()


def _get_all_syn_from_syn_snvs(syn_df: pd.DataFrame, snvs: pd.DataFrame) -> pd.DataFrame:
    return _get_syn_from_snvs(syn_df, snvs, syn_code)


def _get_top_syn_from_mis_snvs(syn_top_df: pd.DataFrame, snvs: pd.DataFrame) -> pd.DataFrame:
    return _get_syn_from_snvs(syn_top_df, snvs, mis_code)


def _get_top_syn_from_non_snvs(syn_top_df: pd.DataFrame, snvs: pd.DataFrame) -> pd.DataFrame:
    return _get_syn_from_snvs(syn_top_df, snvs, non_code)


def _get_snvre_mseq(seq: str, snvre_meta: pd.DataFrame) -> pd.Series:
    def f(item) -> str:
        pos: int = int(item.pos)
        return f"{seq[:pos]}{item.alt}{seq[pos + 3:]}" if pos > 0 else f"{item.alt}{seq[3:]}"

    return snvre_meta[['pos', 'alt']].apply(f, axis=1).astype('string')


def get_fast_snvre(
    syn_df: pd.DataFrame,
    syn_top_df: pd.DataFrame,
    genomic_range: GenomicRange,
    frame: int,
    seq: str,
    snvs: pd.DataFrame
) -> pd.DataFrame:
    """
    Get SNVRE variants as a dataframe given a SNV variant dataframe

    Auxiliary tables syn_df and syn_top_df are computed once per run.

    syn_df: get_all_syn_codons
    syn_top_df: get_top_syn_codons
    snvs: get_fast_exon_snv(..., reset_index=False)
    """

    if not isinstance(snvs.index, pd.CategoricalIndex):
        raise ValueError("Invalid index for SNV metadata table!")

    start_pos: int = genomic_range.start
    end_pos: int = genomic_range.end

    snvs_local: pd.DataFrame = snvs[SNV_SYN_FIELDS].copy()

    # Filter out SNV's that belong to triplets exceeding the target
    inner_cds_start: int = start_pos if frame == 0 else (start_pos + 3 - frame)
    if inner_cds_start != start_pos:
        snvs_local = snvs_local[snvs_local.pos >= np.int32(inner_cds_start)]
    inner_cds_end: int = end_pos - ((end_pos - inner_cds_start + 1) % 3)
    if inner_cds_end != end_pos:
        if inner_cds_end < start_pos:
            return pd.DataFrame([])
        snvs_local = snvs_local[snvs_local.pos <= np.int32(inner_cds_end)]

    # Get synonymous mutations
    syn_meta: pd.DataFrame = _get_all_syn_from_syn_snvs(syn_df, snvs_local)

    # Get top-ranking missense mutation
    mis_meta: pd.DataFrame = _get_top_syn_from_mis_snvs(syn_top_df, snvs_local)

    # Get top-ranking nonsense mutation
    non_meta: pd.DataFrame = _get_top_syn_from_non_snvs(syn_top_df, snvs_local)

    # Merge SNVRE mutations
    meta = pd.concat([
        mut_type_meta
        for mut_type_meta in [syn_meta, mis_meta, non_meta]
        if mut_type_meta.shape[0] != 0
    ], ignore_index=True)

    # TODO: to remove when oligonucleotide names will be generated from metadata directly
    meta.pos -= np.int32(start_pos)

    # Convert parent SNV positions to SNVRE positions
    offset: int = (3 - frame) if frame > 0 else 0
    meta.pos = pd.Series((
        (offset + 3 * ((meta.pos - offset) // 3)) if offset > 0 else
        (3 * (meta.pos // 3))
    ), dtype=np.int32)

    # Drop duplicate alternative codons
    meta = meta.drop_duplicates().reset_index(drop=True)

    # Render triplet replacements
    meta['mseq'] = _get_snvre_mseq(seq, meta)

    # Add variant type
    meta['var_type'] = np.int8(VariantType.SUBSTITUTION.value)

    return meta
