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

import pytest
from valiant.cli_utils import load_codon_table
from valiant.sge_cli import _load_gff_file
from valiant.models.codon_table import CodonTable
from valiant.models.exon import AnnotationRepository
from valiant.models.variant_repository_collection import _load_pam_protection_vcf
from .constants import CODON_TABLE_FP, GTF_SINGLE, GTF_MULTI, PAM_VCF_FP
from .utils import get_data_file_path


def test_load_codon_table_default():

    # Load default codon table
    default_codon_table = load_codon_table(None)
    assert isinstance(default_codon_table, CodonTable)
    aas = set(default_codon_table.amino_acid_symbols)
    assert 'STOP' in aas
    assert len(aas) == 21


def test_load_codon_table_invalid():

    # Load invalid codon table
    with pytest.raises(SystemExit):
        load_codon_table(get_data_file_path(GTF_SINGLE))


def test_load_codon_table_valid():

    # Load valid codon table
    codon_table = load_codon_table(get_data_file_path(CODON_TABLE_FP))

    # Check codon table
    assert isinstance(codon_table, CodonTable)


def test_load_gff_file_none():

    # Load no GTF
    assert _load_gff_file(None) is None


def test_load_gff_file_invalid():
    with pytest.raises(SystemExit):

        # Load invalid GTF
        _load_gff_file(get_data_file_path(GTF_MULTI))


def test_load_gff_file_valid():

    # Load valid GTF
    annotation = _load_gff_file(get_data_file_path(GTF_SINGLE))

    # Check annotation repository
    assert isinstance(annotation, AnnotationRepository)
    assert annotation.cds is not None


def test_load_pam_protection_vcf_none():

    # Load no PAM VCF file
    vr = _load_pam_protection_vcf(frozenset(), None)
    assert vr.sgrna_ids == set()


def test_load_pam_protection_vcf_invalid():

    # Load invalid PAM VCF file
    with pytest.raises(SystemExit):
        _load_pam_protection_vcf(
            frozenset(['sgRNA_1']), get_data_file_path(GTF_SINGLE))


def test_load_pam_protection_vcf_valid():
    sgrna_ids = {'sgRNA_1', 'sgRNA_2', 'sgRNA_3', 'sgRNA_4'}
    pam_fp = get_data_file_path(PAM_VCF_FP)

    # Load valid PAM VCF file
    vr = _load_pam_protection_vcf(sgrna_ids, pam_fp)

    # Check PAM variant repository
    assert vr.sgrna_ids == sgrna_ids
    assert len(vr._variants) == 4
    assert len(vr._ranges) == 5
