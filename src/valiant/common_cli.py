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

from functools import wraps
import click
from .constants import DEFAULT_OLIGO_MAX_LENGTH, DEFAULT_OLIGO_MIN_LENGTH


existing_file = click.Path(exists=True, file_okay=True, dir_okay=False)
writable_dir = click.Path(exists=True, file_okay=False, dir_okay=True, writable=True)


# See discussion: https://github.com/pallets/click/issues/108
def common_params(f):
    @click.argument('oligo-info', type=existing_file)
    @click.argument('ref-fasta', type=existing_file)
    @click.argument('output', type=writable_dir)
    @click.argument('species')
    @click.argument('assembly')
    @click.option('--codon-table', type=existing_file, help="Codon table file path")
    @click.option('--adaptor-5', help="5' adaptor sequence")
    @click.option('--adaptor-3', help="3' adaptor sequence")
    @click.option('--log', default='WARNING', help="Logging level")
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
        return f(*args, **kwargs)
    return wrapper
