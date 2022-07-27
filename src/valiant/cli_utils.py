########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
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

"""
Functions that may force the program to exit, meant to be called by the CLI commands.
"""

import logging
import sys
from typing import Optional

from .models.options import Options
from .constants import DEFAULT_CODON_TABLE_FILE_NAME
from .models.codon_table import CodonTable
from .utils import get_data_file_path


def set_logger(log: str) -> None:
    log_level = log.upper()
    if log_level not in logging._nameToLevel:
        logging.critical("Invalid log level '%s'!" % log)
        sys.exit(1)
    logging.basicConfig(level=logging._nameToLevel[log_level])


def load_codon_table(fp: Optional[str]) -> CodonTable:
    if not fp:
        logging.info("Codon table not specified, the default one will be used.")

    codon_table_fp: str = fp or get_data_file_path(
        DEFAULT_CODON_TABLE_FILE_NAME)

    logging.debug("Loading codon table...")
    try:
        return CodonTable.load(codon_table_fp)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical("Failed to load codon table!")
        sys.exit(1)


def log_excluded_oligo_counts(options: Options, short_oligo_n: int, long_oligo_n: int) -> None:

    # Log number of oligonucleotides discarded due to insufficient length
    if short_oligo_n > 0:
        logging.warning(
            "%d oligonucleotides shorter than %d bases were discarded!" %
            (short_oligo_n, options.oligo_min_length))

    # Log number of oligonucleotides discarded due to excessive length
    if long_oligo_n > 0:
        logging.warning(
            "%d oligonucleotides longer than %d bases were discarded!" %
            (long_oligo_n, options.oligo_max_length))
