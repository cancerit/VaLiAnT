########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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

from enum import Enum
from functools import lru_cache
import os
import pathlib
import re
from typing import List, Type, Tuple, FrozenSet, Iterable
import numpy as np
import pandas as pd
from .constants import SRC_TYPES
from .enums import TargetonMutator

dna_complement_tr_table = str.maketrans('ACGT', 'TGCA')
dna_re = re.compile('^[ACGT]+$')

i8_0: np.int8 = np.int8(0)
i8_1: np.int8 = np.int8(1)
i8_2: np.int8 = np.int8(2)

i8_1_2: List[np.int8] = [i8_1, i8_2]
i8_2_1: List[np.int8] = [i8_2, i8_1]


def get_data_file_path(fp: str) -> str:
    return os.path.join(pathlib.Path(__file__).parent.absolute(), 'data', fp)


def is_dna(seq: str) -> bool:
    return dna_re.match(seq) is not None


def is_strand(strand: str) -> bool:
    return strand == '+' or strand == '-'


def validate_strand(strand: str) -> None:
    if not is_strand(strand):
        raise ValueError(f"Invalid strand '{strand}'!")


def get_region(chromosome: str, start: int, end: int) -> str:
    return f"{chromosome}:{start}-{end}"


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(dna_complement_tr_table)


def parse_list(s: str) -> List[str]:
    return [
        item for item in [
            raw.strip()
            for raw in s.split(',')
        ]
        if item
    ]


def get_smallest_int_type(max_value: int) -> Type:
    if max_value < 0:
        raise ValueError("Only non-negative integers supported!")
    return (
        np.int8 if max_value <= np.iinfo(np.int8).max else
        np.int16 if max_value <= np.iinfo(np.int16).max else
        np.int32 if max_value <= np.iinfo(np.int32).max else
        np.int64
    )


def get_id_column(rown: int) -> np.ndarray:
    return np.arange(rown, dtype=get_smallest_int_type(rown))


def get_constant_category(s: str, n: int, categories: List[str] = None) -> pd.Categorical:
    if categories and s not in categories:
        raise ValueError(f"Invalid category '{s}'!")
    return pd.Categorical([s], categories=categories or [s]).repeat(n)


@lru_cache()
def get_empty_category_column(categories: Iterable[str], n: int) -> pd.Categorical:
    return pd.Categorical(np.empty(n), categories=sorted(categories))


@lru_cache(maxsize=4)
def get_out_of_frame_offset(cds_ext_length: int) -> int:
    return (3 - cds_ext_length) if cds_ext_length > 0 else 0


def get_frame_complement(frame: pd.Series) -> pd.Series:
    return frame.replace(i8_1_2, i8_2_1)


def get_inner_cds_relative_boundaries(seq_len: int, frame: int) -> Tuple[int, int]:
    inner_cds_start: int = 0 if frame == 0 else (3 - frame)
    inner_cds_end: int = seq_len - (seq_len - inner_cds_start) % 3
    return inner_cds_start, inner_cds_end


def get_var_types(var_types: pd.Series) -> List[int]:
    return [
        int(k)
        for k, v in dict(var_types.value_counts()).items()
        if v > 0
    ]


@lru_cache(maxsize=16)
def parse_mutator(s: str) -> TargetonMutator:
    try:
        return TargetonMutator(s)
    except ValueError:
        raise ValueError(f"Invalid mutator '{s}'!")


def parse_mutators(s: str) -> FrozenSet[TargetonMutator]:
    return frozenset(map(parse_mutator, parse_list(s)))


def get_source_type_column(src_type: str, n: int) -> pd.Series:
    return get_constant_category(src_type, n, categories=SRC_TYPES)


def repr_enum_list(enums: Iterable[Enum]) -> str:
    return ', '.join(str(e.value) for e in enums)
