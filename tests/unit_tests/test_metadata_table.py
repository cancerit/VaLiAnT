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
