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

from contextlib import contextmanager
import os
import tempfile
import pandas as pd
from pyranges import PyRanges
import pytest
from valiant.enums import VariantType, VariantClassification
from valiant.loaders.vcf import get_vcf, load_vcf, load_vcf_manifest, normalise_variant_record
from valiant.models.pam_protection import PamVariant, PamProtectionVariantRepository
from valiant.models.variant import VariantRepository
from .utils import get_data_file_path


PAM_VCF_FILE = 'pam.vcf'

VCF_MANIFEST_FIELDS = [
    'vcf_alias',
    'vcf_id_tag',
    'vcf_path'
]

VCF_HEADER = """\
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

DEFAULT_VCF = """\
##fileformat=VCFv4.3
##INFO=<ID={id_tag},Number=1,Type=String,Description="Variant identifier">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrX	41334252	{id_1}	G	C	.	.	{id_tag}=1
chrX	41337416	{id_2}	C	T	.	.	.
chrX	41339064	.	G	A	.	.	{id_tag}=3
"""

class TemporaryVariantFile:
    def __init__(self, alias, variants, id_tag=None):
        self.alias = alias
        self.id_tag = id_tag
        self.variants = variants

        # Create and open temporary VCF file
        with tempfile.NamedTemporaryFile(
            'w+', prefix=self.alias, suffix='.vcf', delete=False
        ) as fh:

            # Store file name
            self.file = fh.name

            # Write header
            fh.write("##fileformat=VCFv4.3\n")
            if self.id_tag:
                fh.write(
                    f"##INFO=<ID={self.id_tag},"
                    "Number=1,Type=Integer,Description=\"Some identifier\">\n")
            fh.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n")

            # Write variant records
            qual = ft = info = '.'
            for id_, chromosome, pos, ref, alt in self.variants:
                if self.id_tag:
                    info = f"{self.id_tag}={id_}"
                    id_ = '.'
                fh.write(f"{chromosome}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\t{ft}\t{info}\n")

    def delete(self):
        os.unlink(self.file)


@contextmanager
def get_temporary_vcf_manifest(temp_vcf_files):
    with tempfile.NamedTemporaryFile('w+') as manifest_f:
        manifest_rows = [
            (vcf.alias, vcf.id_tag or '', vcf.file)
            for vcf in temp_vcf_files
        ]

        write_vcf_manifest(manifest_f, manifest_rows)
        yield manifest_f


pam_vcf_fp = get_data_file_path(PAM_VCF_FILE)


@pytest.mark.parametrize('chromosome', ['chrX'])
def test_load_vcf(chromosome):
    df = load_vcf(pam_vcf_fp, {chromosome: (1, 100000000)})
    if df.shape[0] > 0:
        chr_values = df.Chromosome.unique()
        assert len(chr_values) == 1
        assert chr_values[0] == chromosome


def test_pam_variant_repository_load():
    sgrna_ids = {'sgRNA_1', 'sgRNA_4'}
    vr = PamProtectionVariantRepository(sgrna_ids)

    assert set(vr._variants.keys()) == sgrna_ids

    vr.load(pam_vcf_fp)

    for sgrna_id in sgrna_ids:
        assert sgrna_id in vr._variants
        variants = vr._variants[sgrna_id]
        for variant in variants:
            assert isinstance(variant, PamVariant)


def write_vcf_manifest(f, rows):
    f.write(','.join(VCF_MANIFEST_FIELDS) + '\n')
    for row in rows:
        f.write(','.join(row) + '\n')
    f.write('\n')
    f.seek(0)


def write_dummy_vcf(f, ids, id_tag):
    f.write(DEFAULT_VCF.format(
        id_1=ids[0],
        id_2=ids[1],
        id_tag=id_tag or 'SOMETAG'))
    f.seek(0)


def test_load_vcf_manifest():
    vcf_alias_id_tag = [
        ('clinvar_1', 'ALLELEID'),
        ('clinvar_2', 'ALLELEID'),
        ('gnomad', '')
    ]

    with tempfile.NamedTemporaryFile('w+') as manifest_f:

        tmp_vcfs = [
            tempfile.NamedTemporaryFile(
                'w+', prefix=vcf_alias, suffix='.vcf', delete=False)
            for vcf_alias, id_tag in vcf_alias_id_tag
        ]

        manifest_rows = [
            (vcf_alias, id_tag, f.name)
            for (vcf_alias, id_tag), f in zip(vcf_alias_id_tag, tmp_vcfs)
        ]

        write_vcf_manifest(manifest_f, manifest_rows)

        for (_, id_tag), f in zip(vcf_alias_id_tag, tmp_vcfs):
            write_dummy_vcf(f, [11, 22], id_tag)
            f.close()

        try:
            df = load_vcf_manifest(manifest_f.name, {'chrX': (40000000, 42000000)})

            clinvar_1_ids = df.vcf_var_id[df.vcf_alias == 'clinvar_1'].to_numpy()
            assert clinvar_1_ids[0] == '1'
            assert clinvar_1_ids[1] is pd.NA
            assert clinvar_1_ids[2] == '3'

            clinvar_2_ids = df.vcf_var_id[df.vcf_alias == 'clinvar_2'].to_numpy()
            assert clinvar_2_ids[0] == '1'
            assert clinvar_2_ids[1] is pd.NA
            assert clinvar_2_ids[2] == '3'

            gnomad_ids = df.vcf_var_id[df.vcf_alias == 'gnomad'].to_numpy()
            assert gnomad_ids[0] == '11'
            assert gnomad_ids[1] == '22'
            assert gnomad_ids[2] is pd.NA

        finally:

            # Delete VCF files
            for f in tmp_vcfs:
                os.unlink(f.name)


CHROMOSOME = 'X'
POS = 1000
REF_SEQ = 'AAACCCGGGTTT'

@pytest.mark.parametrize('pos_vcf,ref_vcf,alt_vcf,exp_pos,exp_ref,exp_alt,exp_type', [
    (1000, 'AAAC', 'A', 1001, 'AAC', None, VariantType.DELETION),
    (1000, 'A', 'ATTT', 1001, None, 'TTT', VariantType.INSERTION),
    (1000, 'A', 'T', 1000, 'A', 'T', VariantType.SUBSTITUTION)
])
def test_normalise_variant_record(pos_vcf, ref_vcf, alt_vcf, exp_pos, exp_ref, exp_alt, exp_type):
    with tempfile.NamedTemporaryFile('w', delete=True) as tmp:

        # Write temporary VCF file
        tmp.write(VCF_HEADER)
        tmp.write(f"{CHROMOSOME}\t{POS}\t.\t{ref_vcf}\t{alt_vcf}\t.\t.\t.\n")
        tmp.flush()

        with get_vcf(tmp.name) as vcf:
            for record in vcf.fetch():

                # Check record
                print(record)
                assert record.ref == ref_vcf
                assert len(record.alts) == 1
                assert record.alts[0] == alt_vcf
                assert record.pos == pos_vcf

                # Normalise record
                nr = normalise_variant_record(lambda r: None, record)
                chromosome, start, end, ref, alt, var_id, var_type_sub, var_class = nr
                assert chromosome == CHROMOSOME
                assert start == exp_pos
                assert ref == exp_ref
                assert alt == exp_alt
                assert var_id is None
                assert var_type_sub == exp_type.value
                assert var_class == VariantClassification.CLASSIFIED.value


def test_variant_repository_load():
    region_variant_ids = {
        (CHROMOSOME, 1000, 1100): {1},
        (CHROMOSOME, 2000, 2100): {3, 5},
        (CHROMOSOME, 3000, 3100): set(),
    }

    region_df = PyRanges(df=pd.DataFrame.from_records([
        (chromosome, start - 1, end)
        for chromosome, start, end in region_variant_ids.keys()
    ], columns=['Chromosome', 'Start', 'End']))

    # id, chromosome, pos, ref, alt
    vcf_1 = TemporaryVariantFile('clinvar_1', [
        (1, CHROMOSOME, 1040, 'A', 'C'),
        (2, CHROMOSOME, 900, 'C', 'G')
    ], id_tag='ALLELEID')

    vcf_2 = TemporaryVariantFile('gnomad', [
        (3, CHROMOSOME, 2060, 'A', 'C'),
        (4, CHROMOSOME, 1600, 'C', 'G'),
        (5, CHROMOSOME, 2020, 'CAA', 'C')
    ])

    temp_vcf_files = [vcf_1, vcf_2]

    with get_temporary_vcf_manifest(temp_vcf_files) as manifest_f:
        try:
            vr = VariantRepository.load(manifest_f.name, region_df)
        finally:
            for vcf in temp_vcf_files:
                vcf.delete()

    for region, variant_ids in region_variant_ids.items():
        if variant_ids:
            assert region in vr._region_variants
            assert variant_ids == set(
                int(vr._variants[id_].vcf_variant_id)
                for id_ in vr._region_variants[region]
            )
        else:
            assert region not in vr._region_variants
