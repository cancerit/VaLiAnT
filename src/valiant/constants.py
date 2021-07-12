########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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

from .enums import TargetonMutator

DNA_ALPHABET = 'ACGT'
DNA_ALPHABET_SET = set(DNA_ALPHABET)

DEFAULT_CODON_TABLE_FILE_NAME = 'default_codon_table.csv'

REVCOMP_OLIGO_NAME_SUFFIX = '_rc'
CUSTOM_MUTATOR = 'custom'

DEFAULT_OLIGO_MAX_LENGTH = 300

METADATA_FIELDS = [
    'oligo_name',
    'species',
    'assembly',
    'gene_id',
    'transcript_id',
    'src_type',
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'revc',
    'ref_seq',
    'pam_seq',
    'vcf_alias',
    'vcf_var_id',
    'mut_position',
    'ref',
    'new',
    'ref_aa',
    'alt_aa',
    'mut_type',
    'mutator',
    'oligo_length',
    'mseq'
]

METADATA_FIELDS_SET = set(METADATA_FIELDS)

SRC_TYPES = [
    'ref',
    'cdna'
]

GENERIC_MUTATORS = frozenset([
    TargetonMutator.SNV,
    TargetonMutator.DEL1,
    TargetonMutator.DEL2_0,
    TargetonMutator.DEL2_1
])

CDS_ONLY_MUTATORS = frozenset([
    TargetonMutator.IN_FRAME,
    TargetonMutator.SNV_RE,
    TargetonMutator.STOP,
    TargetonMutator.ALA,
    TargetonMutator.AA
])

MUTATOR_CATEGORIES = sorted([m.value for m in (GENERIC_MUTATORS | CDS_ONLY_MUTATORS)])
