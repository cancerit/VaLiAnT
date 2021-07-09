#!/bin/sh

SCRIPT="${PWD}/../compare_exp_results.sh"
OUT_DIR="output_a"
EXP_OUT_DIR="${OUT_DIR}_exp"

${SCRIPT} ${EXP_OUT_DIR} ${OUT_DIR}
