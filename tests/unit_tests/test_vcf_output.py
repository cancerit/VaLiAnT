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

from tempfile import NamedTemporaryFile
import pandas as pd
import pytest
from valiant.enums import VariantType
from valiant.loaders.vcf import write_vcf, get_vcf
from valiant.models.refseq_repository import ReferenceSequenceRepository
from valiant.models.variant import _get_vcf_record_params_constructor, _get_insertion_record, _get_deletion_record

CHROMOSOME = 'X'
REF_START = 10
REF_END = 15
REF_SEQ = 'AAAGGG'
PREV_NT = 'C'
OLIGO_NAME = '<OLIGO_NAME>'
MUTATOR = '<MUTATOR>'


def get_metadata_table(
    var_type,
    mut_pos,
    ref_seq,
    start,
    ref,
    alt,
    chromosome,
    mutator,
    oligo_name
):
    meta = pd.DataFrame.from_dict({
        'mut_position': [mut_pos],
        'pam_seq': [ref_seq],
        'ref_seq': [ref_seq],
        'ref_start': [start],
        'ref': [ref],
        'new': [alt],
        'ref_chr': [chromosome],
        'mutator': [mutator],
        'oligo_name': [oligo_name],
        'vcf_var_id': [None],
        'var_type': [var_type.value],
        'pam_codon_mask': [0],
        'vcf_alias': [None]
    })

    meta.var_type = meta.var_type.astype('int8')
    meta.mut_position = meta.mut_position.astype('int32')
    meta.ref_start = meta.ref_start.astype('int32')
    meta.ref_chr = meta.ref_chr.astype('category')
    meta.mutator = meta.mutator.astype('string')
    meta.oligo_name = meta.oligo_name.astype('string')
    meta.ref_seq = meta.ref_seq.astype('string')
    meta.pam_seq = meta.pam_seq.astype('string')
    meta.ref = meta.ref.astype('category')
    meta.new = meta.new.astype('category')
    meta.pam_codon_mask = meta.pam_codon_mask.astype('int8')
    meta.vcf_alias = meta.vcf_alias.astype('string')

    return meta


@pytest.fixture(scope='module')
def ref_repository():
    chromosome = CHROMOSOME
    start = REF_START
    end = REF_END
    ref_seq = REF_SEQ
    prev_nt = PREV_NT

    # Populate reference sequence repository
    ref_repository = ReferenceSequenceRepository()
    ref_repository.register_region(chromosome, start, end)
    ref_repository.register_nucleotide(chromosome, start - 1, prev_nt)
    ref_repository.register_sequence(chromosome, start, end, ref_seq)

    yield ref_repository


@pytest.mark.parametrize('var_type,pos,ref,alt,pos_vcf,ref_vcf,alt_vcf', [
    (VariantType.DELETION, 12, 'A', None, 11, 'AA', 'A'),
    (VariantType.DELETION, 12, 'AG', None, 11, 'AAG', 'A'),
    (VariantType.SUBSTITUTION, 12, 'A', 'C', 12, 'A', 'C'),
    (VariantType.SUBSTITUTION, 12, 'AGG', 'TTT', 12, 'AGG', 'TTT'),
    (VariantType.INSERTION, 12, None, 'T', 11, 'A', 'AT'),
    (VariantType.INSERTION, 12, None, 'TT', 11, 'A', 'ATT')
])
def test_get_record(ref_repository, var_type, pos, ref, alt, pos_vcf, ref_vcf, alt_vcf):
    meta = get_metadata_table(
        var_type,
        pos,
        REF_SEQ,
        REF_START,
        ref,
        alt,
        CHROMOSOME,
        MUTATOR,
        OLIGO_NAME)

    to_record_params_f = _get_vcf_record_params_constructor(True)

    def get_record(x):
        return to_record_params_f(x).get_vcf_record(
            ref_repository, True)

    variant = get_record(meta.iloc[0, :])

    stop = (pos_vcf + len(ref_vcf)) if ref is not None else (pos_vcf + 1)

    assert isinstance(variant, dict)
    assert isinstance(variant['alleles'], tuple)
    assert len(variant['alleles']) == 2
    assert variant['alleles'][0] == ref_vcf
    assert variant['alleles'][1] == alt_vcf
    assert variant['contig'] == CHROMOSOME
    assert variant['start'] == pos_vcf - 1
    assert variant['stop'] == stop - 1
    assert isinstance(variant['info'], dict)
    assert variant['info']['SGE_SRC'] == MUTATOR
    assert variant['info']['SGE_OLIGO'] == OLIGO_NAME


def test_get_insertion_record_pos_one():
    ref_seq = 'CGCCGGG'
    pam_seq = 'CGCCGGG'
    alt = 'AAA'
    start_vcf, end_vcf, ref_vcf, alt_vcf, _ = _get_insertion_record(
        ref_repository, CHROMOSOME, 1, alt, 1, ref_seq, pam_seq)

    assert start_vcf == 1
    assert end_vcf == 2

    with NamedTemporaryFile(
        'w+', suffix='.vcf', delete=False
    ) as vcf_fh:
        write_vcf(vcf_fh.name, [CHROMOSOME], [{
            'alleles': (ref_vcf, alt_vcf),
            'contig': CHROMOSOME,
            'start': start_vcf - 1,
            'stop': end_vcf - 1
        }])
        with get_vcf(vcf_fh.name) as vcf:
            record = vcf.fetch().__next__()

        assert record.contig == CHROMOSOME
        assert record.pos == 1
        assert record.ref == pam_seq[0]
        assert len(record.alts) == 1
        assert record.alts[0] == alt + pam_seq[0]
        assert 'END' not in str(record)


def test_get_deletion_record_pos_one():
    ref_seq = 'CGCTCGGG'
    pam_seq = 'CGCTCGGG'
    ref = 'CGC'
    start_vcf, end_vcf, ref_vcf, alt_vcf, _ = _get_deletion_record(
        ref_repository, CHROMOSOME, 1, ref, 1, ref_seq, pam_seq)

    assert start_vcf == 1
    assert end_vcf == start_vcf + len(ref) + 1

    with NamedTemporaryFile(
        'w+', suffix='.vcf', delete=False
    ) as vcf_fh:
        write_vcf(vcf_fh.name, [CHROMOSOME], [{
            'alleles': (ref_vcf, alt_vcf),
            'contig': CHROMOSOME,
            'start': start_vcf - 1,
            'stop': end_vcf - 1
        }])
        with get_vcf(vcf_fh.name) as vcf:
            record = vcf.fetch().__next__()

        assert record.contig == CHROMOSOME
        assert record.pos == 1
        assert record.ref == 'CGCT'
        assert len(record.alts) == 1
        assert record.alts[0] == 'T'
        assert 'END' not in str(record)
