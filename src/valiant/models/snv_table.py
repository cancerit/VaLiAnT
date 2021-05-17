########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
# 
# The Software remains the property of Genome Research Ltd ("GRL").
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
# 
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
# 
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
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
