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

import logging
from contextlib import contextmanager
from typing import Iterable, Generator

from pysam import FastaFile


class SequenceNotFound(Exception):
    pass


# TODO: use explicit index path instead?
def get_fasta_file(fp: str) -> FastaFile:
    try:
        return FastaFile(fp)
    except IOError as ex:
        logging.critical("Failed to load reference file!")
        raise ex


@contextmanager
def open_fasta(fp: str) -> Generator[FastaFile, None, None]:
    ff = get_fasta_file(fp)
    try:
        yield ff
    finally:
        ff.close()


def load_from_multi_fasta(fp: str, ids: Iterable[str]) -> dict[str, str]:
    with open_fasta(fp) as ff:
        try:
            return {
                sid: ff.fetch(sid)
                for sid in ids
            }
        except KeyError as ex:
            logging.critical("Error while loading multi-FASTA file '%s': %s!" % (fp, ex.args[0]))
            raise SequenceNotFound()
