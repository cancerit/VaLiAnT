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

from __future__ import annotations
from dataclasses import dataclass
import logging
from typing import Dict, FrozenSet, Optional
import pandas as pd
from ..snv import _build_base_snv_table, build_snv_table, get_fast_exon_snv, get_all_syn_codons, get_fast_snvre, get_top_syn_codons, build_aa_sub_table, get_fast_aa_subs
from .base import StrandedPositionRange
from .codon_table import CodonTable, STOP_CODE


@dataclass
class AllAminoAcidsTable:
    __slots__ = {'_strand_table'}

    _strand_table: Dict[str, pd.DataFrame]

    def __init__(self, codon_table: CodonTable, strands: FrozenSet[str]) -> None:
        self._strand_table = {
            strand: build_aa_sub_table(codon_table, strand)
            for strand in strands
        }

    def get_subs(self, genomic_range: StrandedPositionRange, frame: int, seq: str) -> pd.DataFrame:
        return get_fast_aa_subs(
            self._strand_table[genomic_range.strand], frame, seq)


@dataclass
class SnvTable:
    __slots__ = {'_codon_table', '_strand_snv_table'}

    _codon_table: CodonTable
    _strand_snv_table: Dict[str, pd.DataFrame]

    def __init__(self, codon_table: CodonTable, strands: FrozenSet[str]) -> None:
        self._codon_table = codon_table

        base_snv_table: pd.DataFrame = _build_base_snv_table()
        self._strand_snv_table = {
            strand: build_snv_table(base_snv_table, codon_table, strand)
            for strand in strands
        }

    @property
    def codon_table(self) -> CodonTable:
        return self._codon_table

    @property
    def strands(self) -> FrozenSet[str]:
        return frozenset(self._strand_snv_table.keys())

    def get_snvs(
        self,
        strand: str,
        cds_seq: str,
        genomic_range: StrandedPositionRange,
        cds_prefix_length: int,
        cds_suffix_length: int,
        reset_index: bool = True
    ) -> pd.DataFrame:
        return get_fast_exon_snv(
            self._strand_snv_table[strand],
            cds_seq,
            genomic_range,
            cds_prefix_length,
            cds_suffix_length,
            reset_index=reset_index)


@dataclass
class SnvReTable:
    __slots__ = {'_snv_table', '_strand_all_syn_table', '_strand_top_syn_table'}

    _snv_table: SnvTable
    _strand_all_syn_table: Dict[str, pd.DataFrame]
    _strand_top_syn_table: Dict[str, pd.DataFrame]

    def __init__(self, snv_table: SnvTable) -> None:
        self._snv_table = snv_table

        # Compute metadata for synonymous mutations
        strand_all_syn_tables: Dict[str, pd.DataFrame] = {
            strand: get_all_syn_codons(self._snv_table.codon_table, strand)
            for strand in snv_table.strands
        }

        # Store metadata for all synonymous mutations
        self._strand_all_syn_table = {
            strand: all_syn_table.drop('codon_rank', axis=1)[all_syn_table.alt_aa != STOP_CODE]
            for strand, all_syn_table in strand_all_syn_tables.items()
        }

        # Filter metadata for synonymous mutations to top ranking codons only
        self._strand_top_syn_table = {
            strand: get_top_syn_codons(all_syn_table)
            for strand, all_syn_table in strand_all_syn_tables.items()
        }

    def get_snvres(self, genomic_range: StrandedPositionRange, frame: int, seq: str, snvs: pd.DataFrame) -> pd.DataFrame:
        strand: str = genomic_range.strand
        return get_fast_snvre(
            self._strand_all_syn_table[strand],
            self._strand_top_syn_table[strand],
            genomic_range,
            frame,
            seq,
            snvs)


@dataclass
class AuxiliaryTables:
    __slots__ = {'codon_table', '_snv_table', '_snvre_table', '_all_aa_table'}

    codon_table: CodonTable
    _snv_table: Optional[SnvTable]
    _snvre_table: Optional[SnvReTable]
    _all_aa_table: Optional[AllAminoAcidsTable]

    def __init__(
        self,
        codon_table: CodonTable,
        strands: FrozenSet[str],
        snv: bool,
        snvre: bool,
        all_aa: bool
    ) -> None:
        self.codon_table = codon_table
        if snv or snvre:
            logging.debug("Initialising SNV tables...")
            self._snv_table = SnvTable(codon_table, strands)
            if snvre:
                logging.debug("Initialising SNVRE tables...")
                self._snvre_table = SnvReTable(self.snv_table)
        if all_aa:
            logging.debug("Initialising all amino acid top ranking tables...")
            self._all_aa_table = AllAminoAcidsTable(codon_table, strands)

    @property
    def snv_table(self) -> SnvTable:
        if not self._snv_table:
            raise RuntimeError("SNV tables not initialised!")
        return self._snv_table

    @property
    def snvre_table(self) -> SnvReTable:
        if not self._snvre_table:
            raise RuntimeError("SNVRE tables not initialised!")
        return self._snvre_table

    @property
    def all_aa_table(self) -> AllAminoAcidsTable:
        if not self._all_aa_table:
            raise RuntimeError("All amino acid top ranking tables not initialised!")
        return self._all_aa_table
