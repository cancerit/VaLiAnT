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

import logging
from typing import Dict, Iterable
from pysam import FastaFile
from ..errors import SequenceNotFound


# TODO: use explicit index path instead?
def get_fasta_file(fp: str) -> FastaFile:
    try:
        return FastaFile(fp)
    except IOError as ex:
        logging.critical("Failed to load reference file!")
        raise ex


def load_from_multi_fasta(fp: str, ids: Iterable[str]) -> Dict[str, str]:
    ff: FastaFile = get_fasta_file(fp)
    try:
        return {
            sid: ff.fetch(sid)
            for sid in ids
        }
    except KeyError as ex:
        logging.critical("Error while loading multi-FASTA file '%s': %s!" % (fp, ex.args[0]))
        raise SequenceNotFound()
    finally:
        ff.close()
