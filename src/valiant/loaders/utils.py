########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023, 2024 Genome Research Ltd
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

from __future__ import annotations

import logging
from enum import IntEnum
from typing import Type

from charset_normalizer import detect

from valiant.uint_range import UIntRange


def detect_encoding(fp: str):
    with open(fp, 'rb') as rfh:
        encoding = detect(rfh.read(10000))['encoding']
    logging.debug("File '%s' encoding: %s." % (fp, encoding))
    return encoding


def parse_list(s: str, delimiter: str = ',', n: int | None = None) -> list[str]:
    ls = [
        item for item in [
            raw.strip()
            for raw in s.split(delimiter)
        ]
        if item
    ]
    if n is not None and len(ls) != n:
        raise ValueError("Unexpected list length!")
    return ls


def get_int_enum(name: str, fields: list[str]) -> Type[IntEnum]:
    return IntEnum(name, [f.upper() for f in fields], start=0)


def parse_uint_range(start: str, end: str) -> UIntRange:
    return UIntRange(int(start), int(end))


def parse_uint_range_from_list(a: list[str], start_field: int, end_field: int) -> UIntRange:
    return parse_uint_range(a[start_field], a[end_field])
