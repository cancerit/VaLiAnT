#!/bin/bash

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
