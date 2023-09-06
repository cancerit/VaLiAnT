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

import pytest
from valiant.enums import TargetonMutator
from valiant.models.base import PositionRange
from valiant.models.refseq_ranges import ReferenceSequenceRanges, ReferenceSequenceRangeCollection, TargetonConfig


def test_reference_sequence_ranges_init():
    r2_mutators = frozenset({TargetonMutator.DEL1})
    mutators = frozenset(), r2_mutators, frozenset()
    chromosome = 'X'
    strand = '+'
    ref_start = 1
    ref_end = 1000
    r2_start = 100
    r2_end = 200
    c1_start = 1
    c1_end = 99
    c2_start = 221
    c2_end = 1000
    r3_start = 201
    r3_end = 220
    sgrna_ids = frozenset({'ID1', 'ID2'})

    # Initialise object
    rsr = ReferenceSequenceRanges.from_config(
        TargetonConfig(
            chromosome,
            strand,
            ref_start,
            ref_end,
            r2_start,
            r2_end,
            (0, 20),
            mutators,
            sgrna_ids))

    # Check reference range
    assert rsr.ref_range.chromosome == chromosome
    assert rsr.ref_range.strand == strand
    assert rsr.ref_range.start == ref_start
    assert rsr.ref_range.end == ref_end

    # Check sgRNA ID's
    assert rsr.sgrna_ids == sgrna_ids

    # Check constant regions
    assert len(rsr._const_regions) == 2
    c1, c2 = rsr._const_regions
    assert c1.chromosome == chromosome
    assert c1.strand == strand
    assert c1.start == c1_start
    assert c1.end == c1_end
    assert c2.chromosome == chromosome
    assert c2.strand == strand
    assert c2.start == c2_start
    assert c2.end == c2_end

    # Check target regions
    assert len(rsr._target_regions) == 3
    r1, r2, r3 = rsr._target_regions
    assert r1 is None
    assert r2.mutators == r2_mutators
    assert r2.genomic_range.chromosome == chromosome
    assert r2.genomic_range.strand == strand
    assert r2.genomic_range.start == r2_start
    assert r2.genomic_range.end == r2_end
    assert not r3.mutators
    assert r3.genomic_range.chromosome == chromosome
    assert r3.genomic_range.strand == strand
    assert r3.genomic_range.start == r3_start
    assert r3.genomic_range.end == r3_end


@pytest.mark.parametrize('s,exp_mutators', [
    ('(1del), (1del,snv), ()', [
        {TargetonMutator.DEL1},
        {TargetonMutator.DEL1, TargetonMutator.SNV},
        set()
    ])
])
def test_reference_sequence_ranges_parse_mutators(s, exp_mutators):
    mutators = ReferenceSequenceRanges.parse_mutator_tuples(s)
    assert mutators == exp_mutators


def test_reference_sequence_range_collection_init():
    r2_mutators = frozenset({TargetonMutator.DEL1})
    mutators = frozenset(), r2_mutators, frozenset()
    chromosome = 'X'
    strand = '+'
    ref_start = 1
    ref_end = 1000
    r2_start = 100
    r2_end = 200
    sgrna_ids = frozenset({'ID1', 'ID2'})

    rsr = ReferenceSequenceRanges.from_config(
        TargetonConfig(
            chromosome,
            strand,
            ref_start,
            ref_end,
            r2_start,
            r2_end,
            (0, 20),
            mutators,
            sgrna_ids))

    rsrc = ReferenceSequenceRangeCollection([rsr])

    # Check reference sequence ranges
    assert len(rsrc._rsrs) == 1
    assert rsrc._rsrs[0] == rsr

    # Check sgRNA ID's
    assert rsrc.sgrna_ids == sgrna_ids

    # Check range tables
    assert len(rsrc._ref_ranges) == 1
    assert len(rsrc._region_ranges) == 4


def get_reference_sequence_ranges(ref_start, ref_end, r2_start, r2_end):
    return ReferenceSequenceRanges.from_config(
        TargetonConfig(
            'X',
            '+',
            ref_start,
            ref_end,
            r2_start,
            r2_end,
            (0, 0),
            (frozenset(), frozenset(), frozenset()),
            frozenset()))


# Region 2 and both constant regions
rsr_c1r2c2 = get_reference_sequence_ranges(100, 400, 200, 299)


@pytest.mark.parametrize('gr,exp', [
    (PositionRange(start, end), exp)
    for start, end, exp in [
        (105, 190, True),  # C1
        (105, 250, True),  # C1 + R2
        (105, 400, True),  # C1 + R2 + C2
        (305, 350, True),  # C2,
        (250, 350, True),  # R2 + C2
        (250, 260, False)  # R2
    ]
])
def test_reference_sequence_ranges_is_range_in_constant_region(gr, exp):

    # Check for overlap with either constant region
    assert rsr_c1r2c2.is_range_in_constant_region(gr) is exp
