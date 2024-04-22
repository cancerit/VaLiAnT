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

import csv

from .strings.codon import Codon
from .codon_table_row import CodonTableRow
from .strings.translation_symbol import TranslationSymbol

RANK_PREFIX = 'RANK'
RANK_PREFIX_LENGTH = len(RANK_PREFIX)
RANK_TOP = 'UT'


def _parse_freq(freq: str) -> float:
    try:
        freq_f = float(freq)
    except ValueError:
        raise ValueError(f"invalid codon frequency '{freq}'")

    if not (0.0 <= freq_f <= 1.0):
        raise ValueError(f"invalid codon frequency '{freq}'")

    return freq_f


def _parse_codon(codon: str) -> Codon:
    try:
        return Codon(codon)
    except ValueError:
        raise ValueError(f"invalid codon '{codon}'")


def _parse_aa(aa: str) -> TranslationSymbol:
    try:
        return TranslationSymbol(aa)
    except ValueError:
        raise ValueError(f"invalid amino acid code '{aa}'")


def _parse_rank(rank: str) -> int:
    if not (rank.startswith(RANK_PREFIX) and len(rank) > RANK_PREFIX_LENGTH):
        raise ValueError(f"invalid codon rank '{rank}'")

    rank_suffix: str = rank[RANK_PREFIX_LENGTH:]
    try:
        return 1 if rank_suffix in RANK_TOP else int(rank_suffix)
    except ValueError:
        raise ValueError(f"invalid codon rank '{rank}'")


def _parse_codon_table_row(it) -> CodonTableRow:
    try:
        codon, aa, freq, rank = it
    except ValueError:
        raise ValueError("invalid number of columns")

    _ = _parse_freq(freq)
    return CodonTableRow(_parse_codon(codon), _parse_aa(aa), _parse_rank(rank))


def load_codon_table_rows(fp: str) -> list[CodonTableRow]:
    try:
        with open(fp) as fh:
            return list(map(_parse_codon_table_row, csv.reader(fh)))
    except ValueError as ex:
        raise ValueError(f"Invalid codon table format: {ex.args[0]}!")
