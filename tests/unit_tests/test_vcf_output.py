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

from tempfile import NamedTemporaryFile
import pandas as pd
import pytest
from valiant.enums import VariantType
from valiant.loaders.vcf import write_vcf, get_vcf
from valiant.models.refseq_repository import ReferenceSequenceRepository
from valiant.models.variant import get_record, _get_insertion_record, _get_deletion_record

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
        'var_type': [var_type.value]
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

    variant = get_record(
        ref_repository, meta.itertuples(index=False).__next__())

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
