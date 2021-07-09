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

EXP_OUT_DIR=$1
OUT_DIR=$2


if command -v md5sum; then
    MD5_CMD='md5sum'
elif command -v md5; then
    MD5_CMD='md5 -r'
else
    echo "Neither md5sum or md5 command found!"
    exit 1
fi

function get_vcfs_md5 () {
    grep -vh '^#' ${1}/*.vcf | sort -s | ${MD5_CMD}
}

function get_csvs_md5 () {
    cat ${1} | sort -s | ${MD5_CMD}
}

function get_unique_csvs_md5 () {
    for f in ${1}/*_unique.csv; do tail -n +2 "$f"; done | sort -s | ${MD5_CMD}
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

get_and_compare_md5s 'get_unique_csvs_md5' "Comparing unique oligonucleotide sequence CSV files" ${EXP_OUT_DIR} ${OUT_DIR}
get_and_compare_md5s 'get_meta_csvs_md5' "Comparing oligonucleotide metadata CSV files" ${EXP_OUT_DIR} ${OUT_DIR}


VCF_FILES=("${EXP_OUT_DIR}"/*.vcf)
if [ -f "${VCF_FILES[0]-}" ]; then
    get_and_compare_md5s 'get_vcfs_md5' "Comparing VCF files" ${EXP_OUT_DIR} ${OUT_DIR}
else
    echo "No VCF files."
fi

META_EXC_FILES=("${EXP_OUT_DIR}"/*_meta_excluded.csv)
if [ -f "${META_EXC_FILES[0]-}" ]; then
    get_and_compare_md5s 'get_meta_excl_csvs_md5' "Comparing excluded oligonucleotide metadata CSV files" ${EXP_OUT_DIR} ${OUT_DIR}
else
    echo "No excluded oligonucleotide metadata files."
fi
