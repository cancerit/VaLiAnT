#!/bin/sh

SCRIPT="${PWD}/../compare_exp_results.sh"
OUT_DIR="brca1_nuc_output"
EXP_OUT_DIR="${OUT_DIR}_exp"

${SCRIPT} ${EXP_OUT_DIR} ${OUT_DIR}
