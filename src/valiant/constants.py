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

from .enums import TargetonMutator

DNA_ALPHABET = 'ACGT'
DNA_ALPHABET_SET = set(DNA_ALPHABET)

DEFAULT_CODON_TABLE_FILE_NAME = 'default_codon_table.csv'

REVCOMP_OLIGO_NAME_SUFFIX = '_rc'
CUSTOM_MUTATOR = 'custom'

MUTATION_TYPE_NON_CDS = 'ncd'

DEFAULT_OLIGO_MAX_LENGTH = 300

# Valid metadata table field names
META_OLIGO_NAME = 'oligo_name'
META_SPECIES = 'species'
META_ASSEMBLY = 'assembly'
META_GENE_ID = 'gene_id'
META_TRANSCRIPT_ID = 'transcript_id'
META_SRC_TYPE = 'src_type'
META_REF_CHR = 'ref_chr'
META_REF_STRAND = 'ref_strand'
META_REF_START = 'ref_start'
META_REF_END = 'ref_end'
META_REVC = 'revc'
META_REF_SEQ = 'ref_seq'
META_PAM_SEQ = 'pam_seq'
META_VCF_ALIAS = 'vcf_alias'
META_VCF_VAR_ID = 'vcf_var_id'
META_MUT_POSITION = 'mut_position'
META_REF = 'ref'
META_NEW = 'new'
META_REF_AA = 'ref_aa'
META_ALT_AA = 'alt_aa'
META_MUT_TYPE = 'mut_type'
META_MUTATOR = 'mutator'
META_OLIGO_LENGTH = 'oligo_length'
META_MSEQ = 'mseq'
META_PAM_MUT_ANNOT = 'pam_mut_annot'
META_PAM_MUT_SGRNA_ID = 'pam_mut_sgrna_id'

# Temporary field names (not included in the final table)
META_VAR_TYPE = 'var_type'

METADATA_FIELDS = [
    META_OLIGO_NAME,
    META_SPECIES,
    META_ASSEMBLY,
    META_GENE_ID,
    META_TRANSCRIPT_ID,
    META_SRC_TYPE,
    META_REF_CHR,
    META_REF_STRAND,
    META_REF_START,
    META_REF_END,
    META_REVC,
    META_REF_SEQ,
    META_PAM_SEQ,
    META_VCF_ALIAS,
    META_VCF_VAR_ID,
    META_MUT_POSITION,
    META_REF,
    META_NEW,
    META_REF_AA,
    META_ALT_AA,
    META_MUT_TYPE,
    META_MUTATOR,
    META_OLIGO_LENGTH,
    META_MSEQ,
    META_PAM_MUT_ANNOT,
    META_PAM_MUT_SGRNA_ID
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
