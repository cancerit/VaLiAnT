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

import csv
from typing import List, Tuple
from ..utils import is_dna


def _parse_codon_table_row(t: Tuple[str, str, str, str]) -> Tuple[str, str, float, int]:
    try:
        codon, aa, freq, rank = t
    except ValueError:
        raise ValueError("invalid number of columns")

    if not (len(codon) == 3 and is_dna(codon)):
        raise ValueError(f"invalid codon code '{codon}'")

    if not (len(aa) == 1 or aa == 'STOP'):
        raise ValueError(f"invalid amino acid code '{aa}'")

    try:
        freq_f = float(freq)
    except ValueError:
        raise ValueError(f"invalid codon frequency '{freq}'")

    if not (0.0 <= freq_f <= 1.0):
        raise ValueError(f"invalid codon frequency '{freq}'")

    RANK_PREFIX = 'RANK'
    rank_prefix_len = len(RANK_PREFIX)

    if not (rank.startswith(RANK_PREFIX) and len(rank) > rank_prefix_len):
        raise ValueError(f"invalid codon rank '{rank}'")

    rank_suffix: str = rank[rank_prefix_len:]
    try:
        rank_i: int = 1 if rank_suffix in 'UT' else int(rank_suffix)
    except ValueError:
        raise ValueError(f"invalid codon rank '{rank}'")

    return codon, aa, freq_f, rank_i


def load_codon_table(fp: str) -> List[Tuple[str, str, float, int]]:

    try:
        with open(fp) as fh:
            return list(map(_parse_codon_table_row, csv.reader(fh)))  # type: ignore
    except ValueError as ex:
        raise ValueError(f"Invalid codon table format: {ex.args[0]}!")
