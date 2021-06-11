#!/bin/sh

EXP_OUT_DIR=$1
OUT_DIR=$2

function get_vcfs_md5 () {
    grep -vh '^#' ${1}/*.vcf | sort | md5sum
}

function get_csvs_md5 () {
    cat ${1} | sort -s | md5sum
}

function get_unique_csvs_md5 () {
    for f in ${1}/*_unique.csv; do tail -n +1 "$f"; done | sort -s | md5sum
}

function get_meta_csvs_md5 () {
    get_csvs_md5 "${1}/*_meta.csv"
}

function get_meta_excl_csvs_md5 () {
    get_csvs_md5 "${1}/*_meta_excluded.csv"
}

function get_and_compare_md5s () {
    local FUNC_NAME="$1"
    local MSG="$2"
    local DIR_A="$3"
    local DIR_B="$4"

    printf '%s... ' "${MSG}"

    MD5=`${FUNC_NAME} "${DIR_A}"`
    EXP_MD5=`${FUNC_NAME} "${DIR_B}"`

    if [ "${MD5}" = "${EXP_MD5}" ]; then
        echo "valid."
    else
        echo "invalid!"
    fi;
}

echo "Comparing results against expected."

get_and_compare_md5s 'get_vcfs_md5' "Comparing VCF files" ${EXP_OUT_DIR} ${OUT_DIR}
get_and_compare_md5s 'get_unique_csvs_md5' "Comparing unique oligonucleotide sequence CSV files" ${EXP_OUT_DIR} ${OUT_DIR}
get_and_compare_md5s 'get_meta_csvs_md5' "Comparing oligonucleotide metadata CSV files" ${EXP_OUT_DIR} ${OUT_DIR}
get_and_compare_md5s 'get_meta_excl_csvs_md5' "Comparing excluded oligonucleotide metadata CSV files" ${EXP_OUT_DIR} ${OUT_DIR}
