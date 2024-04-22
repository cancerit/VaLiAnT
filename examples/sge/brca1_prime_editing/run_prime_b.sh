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

REF=../ref
INPUT=input
TARGETONS=${INPUT}/input_b.txt
REFERENCE_FILE=${REF}/chr17.fa
GTF_FILE=${INPUT}/uniform/ENST00000357654.9.txt
VCF_FILE=${INPUT}/uniform/brca1_vcf_manifest.csv
PAM_FILE=${INPUT}/brca1_prime_ex2_pam.vcf
OUTPUT=output_b
SPECIES=homo_sapiens
ASSEMBLY=GRCh38

# Create output directory
mkdir -p ${OUTPUT}

valiant sge \
    "${TARGETONS}" \
    "${REFERENCE_FILE}" \
    "${OUTPUT}" \
    "${SPECIES}" \
    "${ASSEMBLY}" \
    --gff "${GTF_FILE}" \
    --vcf "${VCF_FILE}" \
    --pam "${PAM_FILE}" \
    --log info \
    --max-length 250 \
    --revcomp-minus-strand \
    --adaptor-5 'AATGATACGGCGACCACCGACACCGAAACTTACTTGCAAAATATGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC' \
    --adaptor-3 'ATTTTGCAAGTAATTTTTCGTATGCCGTCTTCTGCTTG'
