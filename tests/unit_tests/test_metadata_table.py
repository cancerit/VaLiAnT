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

import os
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
from valiant.constants import METADATA_FIELDS
from valiant.models.metadata_table import MetadataTable
from valiant.models.oligo_generation_info import OligoGenerationInfo


METADATA_HEADER = ','.join(METADATA_FIELDS)
UNIQUE_HEADER = 'oligo_name,mseq'

max_oligo_length = 7
partial_meta = pd.DataFrame.from_records([
    (2, 1, 'A', 'C', 'CTGGATTT', 'M', 'L', 'mis', 'snv', 'ATGGATTT', 114, 121, 0, 8, 'human', 'hg38'),
    (2, 1, 'A', 'G', 'GTGGATTT', 'M', 'V', 'mis', 'snv', 'ATGGATTT', 114, 121, 0, 8, 'human', 'hg38'),
    (2, 1, 'A', 'T', 'TTGGATTT', 'M', 'L', 'mis', 'snv', 'ATGGATTT', 114, 121, 0, 8, 'human', 'hg38'),
    (1, 1, 'A', None, 'TGGATTT', None, None, None, '1del', 'ATGGATTT', 114, 121, 0, 7, 'human', 'hg38')
], columns=[
    'var_type',
    'mut_position',
    'ref',
    'new',
    'mseq',
    'ref_aa',
    'alt_aa',
    'mut_type',
    'mutator',
    'ref_seq',
    'ref_start',
    'ref_end',
    'revc',
    'oligo_length',
    'species',
    'assembly'
])

species = 'human'
assembly = 'hg38'


def _check_csv(out_dir, header, rown, fn):
    with open(os.path.join(out_dir, fn)) as fh:
        assert fh.readline().rstrip() == header
        assert len(fh.readlines()) == rown


def test_metadata_table_init():
    mt = MetadataTable(partial_meta.copy(), max_oligo_length)
    assert mt.short_oligo_n == 1
    assert mt.long_oligo_n == 3
    assert np.array_equal(mt.oligo_length_mask.values, [False] * 3 + [True])


def test_metadata_table_from_partial():
    mt = MetadataTable.from_partial(
        species, assembly, partial_meta.copy(), max_oligo_length)
    assert 'species' in mt.metadata
    assert np.array_equal(mt.metadata.species.cat.categories, [species])
    assert 'assembly' in mt.metadata
    assert np.array_equal(mt.metadata.assembly.cat.categories, [assembly])


def test_write_common_files():
    base_fn = 'test_file'
    mt = MetadataTable.from_partial(
        species, assembly, partial_meta.copy(), max_oligo_length)

    with TemporaryDirectory() as out_dir:

        # Write files
        mt.write_common_files(out_dir, base_fn)

        # Check headers
        _check_csv(out_dir, METADATA_HEADER, mt.short_oligo_n, base_fn + '_meta.csv')
        _check_csv(out_dir, METADATA_HEADER, mt.long_oligo_n, base_fn + '_meta_excluded.csv')
        _check_csv(out_dir, UNIQUE_HEADER, mt.short_oligo_n, base_fn + '_unique.csv')


def test_get_info():
    mt = MetadataTable.from_partial(
        species, assembly, partial_meta.copy(), max_oligo_length)
    info = mt.get_info()

    assert isinstance(info, OligoGenerationInfo)
    assert info.long_oligo_n == mt.long_oligo_n
