########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

from valiant.models.alt_seq_builder import AltSeqBuilder
from valiant.models.base import GenomicPosition, GenomicRange
from valiant.models.cds_alt_seq_builder import CdsAltSeqBuilder
from valiant.models.dna_str import DnaStr
from valiant.models.new_pam import PamBgAltSeqBuilder, CdsPamBgAltSeqBuilder
from valiant.models.variant import SubstitutionVariant
from valiant.models.variant_group import VariantGroup

from .utils import get_aux_tables

aux = get_aux_tables()
codon_table = aux.codon_table


def test_alt_seq_builder_get_alt_sub():
    gr = GenomicRange('X', 100, 150, '+')
    i = 5
    pos_a = GenomicPosition(gr.chromosome, gr.start + i)
    pos_b = GenomicPosition(gr.chromosome, gr.start + i + 1)
    ref_nt = 'A'
    var_a = SubstitutionVariant(pos_a, ref_nt, 'T')
    var_b = SubstitutionVariant(pos_b, ref_nt, 'C')
    alt_builder = AltSeqBuilder(gr, 'A' * len(gr), [
        VariantGroup.from_variants([var_a]),
        VariantGroup.from_variants([var_b])
    ])

    assert alt_builder.sequence[i] == ref_nt

    alt = alt_builder.get_alt()
    assert alt[i] == var_a.alt
    assert alt[i + 1] == var_b.alt

    PamBgAltSeqBuilder(alt_builder)


def test_cds_alt_seq_builder_variant_group_codon_clash():
    gr = GenomicRange('X', 100, 123, '+')
    i = 0
    pos_a = GenomicPosition(gr.chromosome, gr.start + i)
    pos_b = GenomicPosition(gr.chromosome, gr.start + i + 1)
    ref_nt = 'A'
    var_a = SubstitutionVariant(pos_a, ref_nt, 'T')
    var_b = SubstitutionVariant(pos_b, ref_nt, 'C')
    alt_builder = CdsAltSeqBuilder(gr, 'A' * len(gr), [
        VariantGroup.from_variants([var_a]),
        VariantGroup.from_variants([var_b])
    ], DnaStr.empty(), DnaStr.empty())

    assert alt_builder.variant_group_codon_clash([0, 1])
    assert alt_builder.variant_group_codon_clash([1, 0])

    var_c = SubstitutionVariant(pos_b, ref_nt, ref_nt)
    assert not var_a.is_frame_shifting
    assert alt_builder.is_variant_nonsynonymous(codon_table, var_a)
    assert not alt_builder.is_variant_nonsynonymous(codon_table, var_c)

    pam_bg = CdsPamBgAltSeqBuilder(alt_builder)
    assert len(pam_bg.bg_variants) == 1
    assert len(pam_bg.pam_variants) == 1
    assert pam_bg.pam_bg_codon_clash
