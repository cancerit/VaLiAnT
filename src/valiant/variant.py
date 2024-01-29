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

from dataclasses import dataclass, replace
from typing import Iterable, NoReturn, TypeVar

from .constants import NT_SNVS
from .enums import VariantType
from .strings.dna_str import DnaStr
from .strings.nucleotide import Nucleotide
from .uint_range import UIntRange
from .utils import get_end


# TODO: drop when frozen dataclasses with slots are fixed in CPython
VariantTuple = tuple[int, DnaStr, DnaStr]


def _raise_no_ref_alt() -> NoReturn:
    raise ValueError("Invalid variant: both REF and ALT are null!")


@dataclass(slots=True)
class Variant:
    pos: int
    ref: DnaStr
    alt: DnaStr

    def __str__(self) -> str:
        return (
            f"{self.pos}del{self.ref}" if not self.alt else
            f"{self.pos}ins{self.alt}" if not self.ref else
            f"{self.pos}{self.ref}>{self.alt}"
        )

    def __post_init__(self) -> None:
        if not self.ref and not self.alt:
            _raise_no_ref_alt()

    @property
    def ref_len(self) -> int:
        return len(self.ref)

    @property
    def alt_len(self) -> int:
        return len(self.alt)

    @property
    def alt_ref_delta(self) -> int:
        return self.alt_len - self.ref_len

    @property
    def ref_end(self) -> int:
        return get_end(self.pos, self.ref_len)

    @property
    def ref_range(self) -> UIntRange:
        return UIntRange(self.pos, self.ref_end)

    @property
    def type(self) -> VariantType:
        if self.ref:
            return (
                VariantType.SUBSTITUTION if self.alt else
                VariantType.DELETION
            )
        if self.alt:
            return VariantType.INSERTION
        _raise_no_ref_alt()

    @property
    def is_insertion(self) -> bool:
        return self.type == VariantType.INSERTION

    @classmethod
    def from_tuple(cls, t: VariantTuple) -> Variant:
        return cls(*t)

    @classmethod
    def from_row(cls, pos: int, ref: str | None, alt: str | None):
        return cls(pos, DnaStr.parse(ref), DnaStr.parse(alt))

    def to_row(self) -> tuple[int, str | None, str | None]:
        return self.pos, self.ref.as_nullable(), self.alt.as_nullable()

    @classmethod
    def get_snvs(cls, pos: int, nt: Nucleotide) -> list[Variant]:
        """Create the three SNV's for a given position"""

        return [
            cls(pos, nt, DnaStr(nt_alt))
            for nt_alt in NT_SNVS[nt]
        ]

    @classmethod
    def get_del(cls, start: int, ref: DnaStr) -> Variant:
        """Create a deletion"""

        return cls(start, ref, DnaStr.empty())

    @classmethod
    def get_ins(cls, start: int, alt: DnaStr) -> Variant:
        """Create an insertion"""

        return cls(start, DnaStr.empty(), alt)

    def clone(self, pos: int | None = None):
        return replace(self, pos=pos if pos is not None else self.pos)

    def offset(self, offset: int):
        return self.clone(pos=self.pos + offset)

    def get_oligo_name_frag(self) -> str:
        pos_frag = (
            str(self.pos) if self.ref_len <= 1 else
            f"{self.pos}_{self.ref_end}"
        )

        match self.type:
            case VariantType.DELETION:
                return pos_frag
            case VariantType.SUBSTITUTION:
                return f"{pos_frag}_{self.ref}>{self.alt}"
            case VariantType.INSERTION:
                return f"{pos_frag}_{self.alt}"
            case VariantType.UNKNOWN:
                raise ValueError("Unknown variant type!")


@dataclass(slots=True)
class RegisteredVariant(Variant):
    id: int


VariantT = TypeVar('VariantT', bound=Variant)


def sort_variants(variants: Iterable[VariantT]) -> list[VariantT]:
    return sorted(variants, key=lambda x: x.pos)
