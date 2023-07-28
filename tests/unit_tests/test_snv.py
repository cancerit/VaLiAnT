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

from io import StringIO
from itertools import chain, repeat
import numpy as np
import pandas as pd
import pytest
from valiant import snv
from valiant.enums import MutationType
from valiant.models.base import GenomicRange
from valiant.models.snv_table import SnvTable, SnvReTable
from .constants import STRANDS
from .utils import get_aux_tables, md5sum, load_csv, trim_cds_extensions

SEQ = 'GGGGGG'
SEQ_LENGTH = len(SEQ)
CDS_SEQ = 'AAAAAAAAAAAA'
CDS_SEQ_LENGTH = len(CDS_SEQ)
SNV_TABLE_ROWN = (4 ** 3) * 9
SNV_BASE_TABLE_COLS = ['triplet', 'pos', 'ref', 'alt', 'alt_triplet']
SNV_TABLE_COLS = SNV_BASE_TABLE_COLS + ['ref_aa', 'alt_aa', 'mut_type']
SNVRE_TABLE_COLS = ['pos', 'ref_aa', 'alt_aa', 'mut_type', 'ref', 'alt', 'mseq', 'var_type']
TOP_SYN_TABLE_COLS = ['syn_triplet', 'ref_aa', 'alt_aa']
ALL_SYN_TABLE_COLS = TOP_SYN_TABLE_COLS

aux = get_aux_tables()
codon_table = aux.codon_table
snv_table = aux.snv_table
snvre_table = aux.snvre_table

syn = MutationType.SYNONYMOUS.value
mis = MutationType.MISSENSE.value
non = MutationType.NONSENSE.value

SNV_DATA = {
    'AAA': (
        ['CAA', 'GAA', 'TAA', 'ACA', 'AGA', 'ATA', 'AAC', 'AAG', 'AAT'],
        np.array([mis, mis, non, mis, mis, mis, mis, syn, mis]),
        np.array([mis, syn, mis, mis, mis, mis, mis, mis, mis])
    ),
    'AGT': (
        ['CGT', 'GGT', 'TGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGC', 'AGG'],
        np.array([mis, mis, mis, mis, mis, mis, mis, syn, mis]),
        np.array([syn, syn, syn, mis, mis, mis, mis, mis, mis])
    )
}

STRAND_LABELS = {
    '+': 'plus',
    '-': 'minus'
}


def normalise_snvre_table(df):
    return df.sort_values(['pos', 'alt_aa', 'alt'], ignore_index=True).reset_index(drop=True)


def get_snvre_fn(seq, strand):
    return f"{seq.lower()}_{STRAND_LABELS[strand]}_snvre.csv"


def load_snvre_table(fp):
    df = load_csv(fp)
    amino_acid_symbols = codon_table.amino_acid_symbols

    df.pos = df.pos.astype('int32')
    df.pos = df.pos.astype('int32')
    df.ref_aa = snv.as_amino_acids(df.ref_aa, amino_acid_symbols)
    df.alt_aa = snv.as_amino_acids(df.alt_aa, amino_acid_symbols)
    df.mut_type = df.mut_type.astype('int8')
    df.ref = snv.as_triplets(df.ref)
    df.alt = snv.as_triplets(df.alt)
    df.mseq = df.mseq.astype('string')
    df.var_type = df.var_type.astype('int8')
    return normalise_snvre_table(df)


snvre_ref = {
    (seq, strand): load_snvre_table(get_snvre_fn(seq, strand))
    for strand, seq in [
        ('+', 'AAAAGT'),
        ('-', 'AAAACT'),
        ('-', 'TTTACT')
    ]
}


def cds_seq_to_genomic_range(cds_seq, strand, plen, slen, offset=10):
    start = offset + plen
    end = offset + len(cds_seq) - slen - 1
    return GenomicRange('X', start, end, strand)


def validate_base_snv_table_dtypes(df):
    assert df.triplet.dtype == 'category'
    assert df.pos.dtype == 'int32'
    assert df.ref.dtype == 'category'
    assert df.alt.dtype == 'category'
    assert df.alt_triplet.dtype == 'category'


def validate_base_snv_table(base_snv_table):
    assert list(base_snv_table.columns) == SNV_BASE_TABLE_COLS
    assert isinstance(base_snv_table.index, pd.CategoricalIndex)
    validate_base_snv_table_dtypes(base_snv_table)
    assert base_snv_table.shape[0] == SNV_TABLE_ROWN


def validate_snv_table_dtypes(df):
    validate_base_snv_table_dtypes(df)
    assert df.ref_aa.dtype == 'category'
    assert df.alt_aa.dtype == 'category'
    assert df.mut_type.dtype == 'int8'


def validate_snv_table(df):
    assert list(df.columns) == SNV_TABLE_COLS
    assert isinstance(df.index, pd.CategoricalIndex)
    validate_snv_table_dtypes(df)
    assert df.shape[0] == SNV_TABLE_ROWN


def validate_snv_meta(seq_length, df):
    assert list(df.columns) == SNV_TABLE_COLS
    validate_snv_table_dtypes(df)
    assert df.shape[0] == 3 * seq_length


def validate_snvre_meta(df):
    assert list(df.columns) == SNVRE_TABLE_COLS
    assert isinstance(df.index, pd.RangeIndex)
    assert df.pos.dtype == 'int32'
    assert df.ref_aa.dtype == 'category'
    assert df.alt_aa.dtype == 'category'
    assert df.mut_type.dtype == 'int8'
    assert df.ref.dtype == 'category'
    assert df.alt.dtype == 'category'
    assert df.mseq.dtype == 'string'
    assert df.var_type.dtype == 'int8'


def test_build_base_snv_table():
    base_snv_table = snv._build_base_snv_table()
    validate_base_snv_table(base_snv_table)


def validate_syn_table_dtypes(df):
    assert df.syn_triplet.dtype == 'category'
    assert df.ref_aa.dtype == 'category'
    assert df.alt_aa.dtype == 'category'


def validate_all_syn_table(df):
    assert list(df.columns) == ALL_SYN_TABLE_COLS
    validate_syn_table_dtypes(df)
    assert df.shape[0] == 174


def validate_top_syn_table(df):
    assert list(df.columns) == TOP_SYN_TABLE_COLS
    validate_syn_table_dtypes(df)
    assert df.shape[0] == 62


@pytest.mark.parametrize('strand,csv_md5', [
    ('+', '70d3994fbf04597059c40dd849cc6279'),
    ('-', '82156651904a1dd71a54a88c15f7f1f6')
])
def test_build_snv_table(strand, csv_md5):

    # Build SNV table
    base_snv_table = snv._build_base_snv_table()
    snv_table = snv.build_snv_table(base_snv_table, codon_table, strand)

    # Check table format
    validate_snv_table(snv_table)

    # Make the representation deterministic
    snv_table = snv_table.sort_values(['pos', 'ref', 'alt'])

    # Check table contents via CSV representation
    with StringIO() as buffer:
        snv_table.to_csv(buffer)
        assert md5sum(buffer.getvalue()) == csv_md5


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('plen', [0, 1, 2])
@pytest.mark.parametrize('slen', [0, 1, 2])
def test_get_fast_exon_snv(strand, plen, slen):

    # Build SNV table
    snv_table = snv.build_snv_table(
        snv._build_base_snv_table(), codon_table, strand)

    # Generate genomic range
    gr = cds_seq_to_genomic_range(CDS_SEQ, strand, plen, slen)
    target_length = len(gr)

    # Generate SNV metadata table
    snvs = snv.get_fast_exon_snv(
        snv_table, CDS_SEQ, gr, plen, slen, reset_index=False)

    # Check table format
    validate_snv_meta(target_length, snvs)

    # Check table contents
    assert list(snvs.pos) == list(
        chain.from_iterable(
            repeat(x, 3)
            for x in range(gr.start, gr.end + 1)
        ))


@pytest.mark.parametrize('records', [[
    (0, 'AAT'),
    (0, 'AAC'),
    (0, 'AAG'),
    (3, 'CAA')
]])
def test_get_snvre_mseq(records):
    meta = pd.DataFrame.from_records(records, columns=['pos', 'alt'])
    mseqs = snv._get_snvre_mseq(SEQ, meta)
    assert isinstance(mseqs, pd.Series)
    assert mseqs.dtype == 'string'
    for (codon_pos, alt), mseq in zip(records, mseqs):
        assert len(mseq) == SEQ_LENGTH
        assert mseq[codon_pos:codon_pos + 3] == alt


def test_snv_table():

    # Build SNV table
    snv_table = SnvTable(codon_table, STRANDS)

    # Check SNV table attributes
    assert snv_table.codon_table == codon_table
    assert snv_table.strands == STRANDS
    for strand in STRANDS:
        validate_snv_table(snv_table._strand_snv_table[strand])


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('plen', [0, 1, 2])
@pytest.mark.parametrize('slen', [0, 1, 2])
def test_snv_table_get_snvs(strand, plen, slen):

    # Build SNV table
    snv_table = SnvTable(codon_table, STRANDS)

    # Generate genomic range
    gr = cds_seq_to_genomic_range(CDS_SEQ, strand, plen, slen)

    # Generate SNV metadata table
    snv_meta = snv_table.get_snvs(
        strand, CDS_SEQ, gr, plen, slen, reset_index=True)

    # Check table format
    validate_snv_meta(len(gr), snv_meta)


@pytest.mark.parametrize('triplet,alt_triplets,mut_types_plus,mut_types_minus', [
    (triplet, alt_triplets, mut_types_plus, mut_types_minus)
    for triplet, (alt_triplets, mut_types_plus, mut_types_minus) in SNV_DATA.items()
])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_snv_table_get_snvs_triplet(triplet, alt_triplets, strand, mut_types_plus, mut_types_minus):
    gr = GenomicRange('X', 10, 12, strand)
    tr = codon_table.translate if strand == '+' else codon_table.translate_rc
    mut_types = mut_types_plus if strand == '+' else mut_types_minus

    # Retrieve SNV metadata
    snv_meta = snv_table.get_snvs(strand, triplet, gr, 0, 0, reset_index=True)
    snv_meta = snv_meta.sort_values(['pos', 'alt']).reset_index(drop=True)

    # Check amino acid change
    ref_aas = snv_meta.ref_aa.unique().astype('string')
    assert len(ref_aas) == 1
    assert ref_aas[0] == tr(triplet)
    assert np.array_equal(snv_meta.alt_aa, np.array(list(map(tr, alt_triplets))))
    assert np.array_equal(snv_meta.mut_type, mut_types)


def test_snvre_table():

    # Build SNVRE table
    snvre_table = SnvReTable(SnvTable(codon_table, STRANDS))

    # Check SNVRE table attributes
    for t in snvre_table._strand_all_syn_table.values():
        validate_all_syn_table(t)
    for t in snvre_table._strand_top_syn_table.values():
        validate_top_syn_table(t)


@pytest.mark.parametrize('strand,cds_seq,sub_seq', [
    ('+', 'GGGAAAAGTGGG', 'AAAAGT'),
    ('-', 'CCCTTTACTCCC', 'TTTACT'),
    ('-', 'GGGAAAACTGGG', 'AAAACT')
])
@pytest.mark.parametrize('plen', [0, 1, 2])
@pytest.mark.parametrize('slen', [0, 1, 2])
def test_snvre_table_get_snvres(strand, cds_seq, sub_seq, plen, slen):

    cds_seq_length = len(cds_seq)

    # Extract target from CDS sequence
    seq = trim_cds_extensions(cds_seq, plen, slen)
    seq_length = len(seq)

    # Generate genomic range
    gr = cds_seq_to_genomic_range(cds_seq, strand, plen, slen)

    # Validate inputs
    assert cds_seq_length % 3 == 0
    assert seq_length == cds_seq_length - plen - slen
    assert seq_length == len(gr)

    # Locate subsequence
    sub_seq_start = seq.index(sub_seq)
    sub_seq_end = sub_seq_start + len(sub_seq)

    # Generate SNV metadata table
    snv_meta = snvre_table._snv_table.get_snvs(
        strand, cds_seq, gr, plen, slen, reset_index=False)

    # Generate SNVRE metadata table
    snvre_meta = snvre_table.get_snvres(gr, plen, seq, snv_meta)

    # Check table format
    validate_snvre_meta(snvre_meta)

    # Check all positions are within the sequence genomic range
    snvre_meta.pos += np.int32(gr.start)
    assert (
        (snvre_meta.pos >= np.int32(gr.start))
        & (snvre_meta.pos <= np.int32(gr.end))
    ).all()
    snvre_meta.pos -= np.int32(gr.start)

    # Check all mutated sequences have the same length of the target
    assert snvre_meta.mseq.str.len().eq(seq_length).all()

    # Remove variants outside of the reference subsequence range
    snvre_meta = snvre_meta[
        (snvre_meta.pos >= np.int32(sub_seq_start))
        & (snvre_meta.pos < np.int32(sub_seq_end))
    ]

    # Remove subsequence range offset from positions
    snvre_meta.pos -= np.int32(sub_seq_start)

    # Trim mutated sequence to reference subsequence range
    snvre_meta.mseq = snvre_meta.mseq.str.slice(sub_seq_start, sub_seq_end)

    # Compare the result subset with the reference results
    snvre_meta_norm = normalise_snvre_table(snvre_meta)

    assert snvre_meta_norm.equals(snvre_ref[sub_seq, strand])
