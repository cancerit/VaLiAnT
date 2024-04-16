########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2024 Genome Research Ltd
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

from pysam import FastaFile

from .seq import Seq
from .strings.dna_str import DnaStr
from .uint_range import UIntRange


def fetch_sequence(fa: FastaFile, contig: str, r: UIntRange) -> Seq:
    logging.debug("Fetching reference sequence at %s:%d-%d." % (contig, r.start, r.end))

    assert r.start > 0
    seq = DnaStr(fa.fetch(reference=contig, start=r.start - 1, end=r.end).upper())

    return Seq(r.start, seq)
