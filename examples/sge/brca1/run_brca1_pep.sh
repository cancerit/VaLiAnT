#!/bin/bash

########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2024 Genome Research Ltd
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

REF=${PWD}/../ref

# Create output directory
mkdir -p brca1_pep_output

# Generate library
valiant sge \
    parameter_input_files/brca1_pep_targeton_input.txt \
    "${REF}/chr17.fa" \
    brca1_pep_output \
    'homo sapiens' \
    'GRCh38' \
    --pam parameter_input_files/brca1_protection_edits.vcf \
    --revcomp-minus-strand \
    --adaptor-5 AATGATACGGCGACCACCGA \
    --adaptor-3 TCGTATGCCGTCTTCTGCTTG \
    --gff reference_input_files/ENST00000357654.9.gtf
