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

import os
import pathlib
import re
from dataclasses import fields
from itertools import groupby
from typing import Callable, Iterable

from .constants import DATA_PATH, CODON_TABLE_FN, DDL_FN

dna_re = re.compile('^[ACGT]*$')
dna_complement_tr_table = str.maketrans('ACGT', 'TGCA')


def is_dna(s: str) -> bool:
    return dna_re.match(s) is not None


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(dna_complement_tr_table)


def parse_opt_int_group(m: re.Match, i: int) -> int:
    g = m.group(i)
    return int(g) if g else 0


def safe_group_by(a: Iterable, k: Callable):
    return groupby(sorted(a, key=k), key=k)


def get_data_file_path(fp):
    return os.path.join(pathlib.Path(__file__).parent.absolute(), DATA_PATH, fp)


def get_default_codon_table_path() -> str:
    return get_data_file_path(CODON_TABLE_FN)


def get_ddl_path() -> str:
    return get_data_file_path(DDL_FN)


def has_duplicates(items: list) -> bool:
    return len(set(items)) != len(items)


def get_dataclass_fields(cls) -> list[str]:
    return [f.name for f in fields(cls)]


def get_enum_values(cls) -> list:
    return [x.value for x in cls]


def get_not_none(it):
    return [x for x in it if x is not None]


def is_adaptor_valid(adaptor: str | None) -> bool:
    return adaptor is None or is_dna(adaptor)


def get_codon_offset_complement(offset: int) -> int:
    match offset:
        case 0:
            return 0
        case 1:
            return 2
        case 2:
            return 1
        case _:
            raise ValueError("Invalid codon offset!")


def get_cds_ext_3_length(frame: int, length: int) -> int:
    """
    Calculate how many nucleotides are missing from the last codon
    given the reading frame and the length of the sequence
    """

    return (3 - (length + frame) % 3) % 3


def has_duplicates(items: list[int]) -> bool:
    return len(set(items)) != len(items)
