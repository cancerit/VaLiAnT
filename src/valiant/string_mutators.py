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

from typing import List, Tuple, Set


def _delete_non_overlapping(seq: str, count: int, offset: int = 0) -> List[Tuple[int, str, str]]:
    return [
        (i, seq[i:i + count], seq[:i] + seq[i + count:])
        for i in range(offset, count * ((len(seq) - offset) // count), count)
    ]


def delete_1(seq: str) -> List[Tuple[int, str, str]]:
    return [
        (i, seq[i], seq[:i] + seq[i + 1:])
        for i in range(len(seq))
    ]


def delete_non_overlapping_2(seq: str, offset: int) -> List[Tuple[int, str, str]]:
    return _delete_non_overlapping(seq, 2, offset=offset)


def delete_non_overlapping_3(seq: str) -> List[Tuple[int, str, str]]:
    return _delete_non_overlapping(seq, 3)


def delete_non_overlapping_3_offset(seq: str, start_offset: int, end_offset: int) -> List[Tuple[int, str, str]]:
    if start_offset < 0 or end_offset < 0:
        raise ValueError("Invalid offset!")

    if start_offset == 0 and end_offset == 0:
        return delete_non_overlapping_3(seq)

    prefix = seq[:start_offset]

    if end_offset > 0:
        suffix = seq[-end_offset:]
        subseq = seq[start_offset:-end_offset]
    else:
        suffix = ''
        subseq = seq[start_offset:]

    return [
        (i + start_offset, ref, f"{prefix}{mseq}{suffix}")
        for i, ref, mseq in delete_non_overlapping_3(subseq)
    ]


def replace_single_nucleotides(seq: str, alphabet: Set[str]) -> List[Tuple[int, str, str, str]]:
    return [
        (i, seq[i], nt, f"{seq[:i]}{nt}{seq[i + 1:]}")
        for i in range(len(seq))
        for nt in alphabet
        if nt != seq[i]
    ]


def replace_nucleotides(seq: str, offset: int, ref: str, alt: str) -> str:
    n: int = len(seq)
    if offset >= n:
        raise ValueError("Substitution out of range!")
    return f"{seq[:offset]}{alt}{seq[offset + len(ref):]}"


# TODO: raise exception on out-of-range (after excluding in variant filtering)
def delete_nucleotides(seq: str, offset: int, ref: str) -> str:
    n: int = len(seq)
    if offset >= n:
        return seq
    return f"{seq[:offset]}{seq[len(ref) + offset:]}"


def insert_nucleotides(seq: str, offset: int, alt: str) -> str:
    n: int = len(seq)
    if offset > n:
        raise ValueError("Insertion out of range!")
    return f"{seq[:offset]}{alt}{seq[offset:]}" if offset < n else (seq + alt)


def replace_codons_const(seq: str, offset: int, codon: str) -> List[Tuple[int, str, str, str]]:
    return [
        (i, seq[i:i + 3], codon, f"{seq[:i]}{codon}{seq[i + 3:]}")
        for i in range(offset, offset + (3 * ((len(seq) - offset) // 3)), 3)
        if seq[i:i + 3] != codon
    ]
