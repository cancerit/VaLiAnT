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

from array import array


def get_zero_array(t: str, length: int) -> array:
    k = array(t).itemsize
    return array(t, bytes(k * length))


def get_u8_array(n: int) -> array:
    return get_zero_array('B', n)


def get_u32_array(n: int) -> array:
    return get_zero_array('I', n)


def get_prev_index(a: array, i: int, value: int) -> int | None:
    j = i - 1
    while j >= 0:
        if a[j] == value:
            return j
        j -= 1
    return None


def get_next_index(a: array, i: int, value: int) -> int | None:
    try:
        return a.index(value, i + 1)
    except ValueError:
        return None
