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

from .annotation import Annotation
from .seq import Seq
from .strings.dna_str import DnaStr
from .uint_range import UIntRange


def fetch_sequence(fa: FastaFile, contig: str, r: UIntRange) -> Seq:
    logging.debug("Fetching reference sequence at %s:%d-%d." % (contig, r.start, r.end))

    assert r.start > 0
    seq = DnaStr(fa.fetch(reference=contig, start=r.start - 1, end=r.end).upper())

    return Seq(r.start, seq)


def get_background_context_range(
    targeton_ranges: list[UIntRange],
    annot: Annotation | None
) -> UIntRange:
    """Get the minimum range in which to apply background variants"""

    # Identify minimal context required
    start = min(t.start for t in targeton_ranges)
    end = max(t.end for t in targeton_ranges)

    if annot:
        start = min(start, annot.cds_start)
        end = max(end, annot.cds_end)

    # Add the preceding nucleotide (for the purposes of VCF indel formatting)
    if start > 1:
        start -= 1

    assert start != 0
    return UIntRange(start, end)


# def select_ppes_in_range(self, sgrna_ids: frozenset[str]) -> GetVariantsInRangeCallable:

#     # Build query (add sgrna ID filter to positional filter)
#     query = SqlQuery(sql_and([
#         sql_select_ppes_in_range,
#         sql_eq_or_in_str_list('sgrna_id', list(sgrna_ids))
#     ]))

#     return VariantSelectStart(query).select_in_range
