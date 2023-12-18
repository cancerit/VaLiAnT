########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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

from functools import wraps
import logging
import os
import sys
import click

from .constants import DEFAULT_OLIGO_MAX_LENGTH, DEFAULT_OLIGO_MIN_LENGTH, OUTPUT_CONFIG_FILE_NAME
from .config import BaseConfig
from .errors import InvalidConfig
from .oligo_generation_info import OligoGenerationInfo
from .main_config import get_main_config_from_config


existing_file = click.Path(exists=True, file_okay=True, dir_okay=False)
writable_dir = click.Path(exists=True, file_okay=False, dir_okay=True, writable=True)


def set_logger(ctx: click.Context, param: click.Parameter, value: str) -> None:
    logging.basicConfig(level=logging._nameToLevel[value.upper()])


# See discussion: https://github.com/pallets/click/issues/108
def common_params(f):
    @click.argument('oligo_info_fp', type=existing_file, metavar='OLIGO_INFO')
    @click.argument('ref_fasta_fp', type=existing_file, metavar='REF_FASTA')
    @click.argument('output_dir', type=writable_dir, metavar="OUTPUT")
    @click.argument('species')
    @click.argument('assembly')
    @click.option('--codon-table', 'codon_table_fp', type=existing_file, help="Codon table file path")
    @click.option('--adaptor-5', help="5' adaptor sequence")
    @click.option('--adaptor-3', help="3' adaptor sequence")
    @click.option(
        '--log',
        default='WARNING',
        type=click.Choice(list(logging._nameToLevel.keys()), case_sensitive=False),
        callback=set_logger,
        help="Logging level")
    @click.option(
        '--max-length',
        type=int,
        default=DEFAULT_OLIGO_MAX_LENGTH,
        help="Maximum oligonucleotide length")
    @click.option(
        '--min-length',
        type=int,
        default=DEFAULT_OLIGO_MIN_LENGTH,
        help="Minimum oligonucleotide length")
    @wraps(f)
    def wrapper(*args, **kwargs):
        kwargs.pop('log', None)

        try:
            return f(*args, **kwargs)
        except InvalidConfig:
            logging.critical("Invalid configuration!")
            sys.exit(1)

    return wrapper


def log_excluded_oligo_counts(config: BaseConfig, short_oligo_n: int, long_oligo_n: int) -> None:

    # Log number of oligonucleotides discarded due to insufficient length
    if short_oligo_n > 0:
        logging.warning(
            "%d oligonucleotides shorter than %d bases were discarded!" %
            (short_oligo_n, config.min_length))

    # Log number of oligonucleotides discarded due to excessive length
    if long_oligo_n > 0:
        logging.warning(
            "%d oligonucleotides longer than %d bases were discarded!" %
            (long_oligo_n, config.max_length))


def finalise(config: BaseConfig, stats: OligoGenerationInfo) -> None:
    """Common operations to be performed at the end independently of the execution mode"""

    # Log the count of oligonucleotides excluded by the length filter
    log_excluded_oligo_counts(config, stats.short_oligo_n, stats.long_oligo_n)

    config_fp: str = config.get_output_file_path(OUTPUT_CONFIG_FILE_NAME)

    try:

        # Generate full configuration
        main_config = get_main_config_from_config(config)

        # Write configuration to file
        main_config.write(config_fp)

    except InvalidConfig:
        logging.error("Failed to generate configuration!")

    except (PermissionError, IsADirectoryError):
        logging.error("Failed to write configuration to '%s'!" % config_fp)
