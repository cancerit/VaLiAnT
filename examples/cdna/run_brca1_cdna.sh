#!/bin/bash

INPUT=${PWD}/brca1_input_files
TARGETONS=${INPUT}/brca1_cdna_targeton.tsv
FASTA=${INPUT}/brca1_203_et_al.fa
ANNOT=${INPUT}/brca1_cdna_annot.tsv
OUTPUT=${PWD}/brca1_cdna_output
SPECIES='homo sapiens'
ASSEMBLY='GRCh38'

# Create output directory
mkdir -p ${OUTPUT}

ls ${OUTPUT}

valiant cdna \
    "${TARGETONS}" \
    "${FASTA}" \
    "${OUTPUT}" \
    "${SPECIES}" \
    "${ASSEMBLY}" \
    --annot "${ANNOT}" \
    --adaptor-5 'GGGAAA' \
    --adaptor-3 'AAAGGG' \
    --log debug
