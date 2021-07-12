#!/bin/bash

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

INPUT=${PWD}/input
TARGETONS=${INPUT}/cdna_targeton.tsv
FASTA=${INPUT}/BRCA1_NP_009225_1_pCW57_1.fa
ANNOT=${INPUT}/cdna_annot.tsv
OUTPUT=${PWD}/output
SPECIES=human
ASSEMBLY=pCW57.1

# Create output directory
mkdir -p ${OUTPUT}

valiant cdna \
    "${TARGETONS}" \
    "${FASTA}" \
    "${OUTPUT}" \
    "${SPECIES}" \
    "${ASSEMBLY}" \
    --annot "${ANNOT}" \
    --adaptor-5 'AATGATACGGCGACCACCGA' \
    --adaptor-3 'TCGTATGCCGTCTTCTGCTTG' \
