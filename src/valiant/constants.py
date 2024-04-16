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

NTS = set('ACGT')

NT_SNVS = {
    nt: sorted(NTS - {nt})
    for nt in NTS
}

# Stop symbol (codon table)
STOP = 'STOP'

# Path to the package data directory
DATA_PATH = 'data'

# Database DDL script file name
DDL_FN = 'ddl.sql'

# Default codon table file name
CODON_TABLE_FN = 'default_codon_table.csv'

# Output configuration file name
OUTPUT_CONFIG_FILE_NAME = 'config.json'

# Output SGE sequences file name
OUTPUT_REF_SEQUENCES = 'ref_sequences.csv'

# Default parameters
DEFAULT_OLIGO_MAX_LENGTH = 300
DEFAULT_OLIGO_MIN_LENGTH = 1

REVCOMP_OLIGO_NAME_SUFFIX = '_rc'
