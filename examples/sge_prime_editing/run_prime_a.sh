#!/bin/bash

INPUT=${PWD}/input
TARGETONS=${INPUT}/input_a.txt
REFERENCE_FILE=${INPUT}/uniform/chr17.fa
GTF_FILE=${INPUT}/uniform/ENST00000357654.9.txt
VCF_FILE=${INPUT}/uniform/brca1_vcf_manifest.csv
PAM_FILE=${INPUT}/brca1_prime_ex2_pam.vcf
OUTPUT=${PWD}/output_a
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
    --adaptor-5 'AATGATACGGCGACCACCGACACCGTGCTAGTCTGGAGTTGATCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC' \
    --adaptor-3 'TCAACTCCAGACTTTTTTCGTATGCCGTCTTCTGCTTG'