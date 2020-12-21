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
from itertools import groupby
from typing import Dict, List, Tuple
from ..globals import TRIPLET_RCS

STOP_CODE = 'STOP'


class CodonTable:
    __slots__ = {'_codon2aa', '_aa2codons', '_codon2syn', '_codonrc2syn'}

    def __init__(self, codon_table: List[Tuple[str, str, float, int]]) -> None:

        # Codon -> Amino acid
        self._codon2aa: Dict[str, str] = {
            codon: aa
            for codon, aa, _, _ in codon_table
        }

        # Amino acid -> Codons (rank sorted)
        self._aa2codons: Dict[str, List[str]] = {
            aa: [t[0] for t in sorted(ts, key=lambda x: x[3])]
            for aa, ts in groupby(codon_table, key=lambda t: t[1])
        }

        # Codon -> Synonymous codons (rank sorted)
        self._codon2syn: Dict[str, List[str]] = {
            codon: [c for c in self._aa2codons[self._codon2aa[codon]] if c != codon]
            for codon in self._codon2aa.keys()
        }

        # Codon -> Synonymous codons of reverse complement (rank sorted)
        self._codonrc2syn: Dict[str, List[str]] = {
            codon: [
                c
                for c in self._aa2codons[self._codon2aa[TRIPLET_RCS[codon]]]
                if c != TRIPLET_RCS[codon]
            ]
            for codon in self._codon2aa.keys()
        }

    @property
    def amino_acid_symbols(self) -> List[str]:
        return list(self._aa2codons.keys())

    def get_top_ranking_codon_table(self) -> List[Tuple[str, str]]:
        return [(aa, codons[0]) for aa, codons in self._aa2codons.items()]

    def get_top_ranking_codon_table_rc(self) -> List[Tuple[str, str]]:
        return [(aa, TRIPLET_RCS[codons[0]]) for aa, codons in self._aa2codons.items()]

    def get_top_ranking_codon(self, aa: str) -> str:
        return self._aa2codons[aa][0]

    def get_top_ranking_codon_rc(self, aa: str) -> str:
        return TRIPLET_RCS[self._aa2codons[aa][0]]

    def translate(self, codon: str) -> str:
        return self._codon2aa[codon]

    def translate_rc(self, codon: str) -> str:
        return self._codon2aa[TRIPLET_RCS[codon]]

    def get_translation_table(self) -> List[Tuple[str, str]]:
        return list(self._codon2aa.items())

    def get_translation_table_rc(self) -> List[Tuple[str, str]]:
        return [
            (TRIPLET_RCS[codon], aa)
            for codon, aa in self._codon2aa.items()
        ]

    @classmethod
    def load(cls, fp: str) -> CodonTable:
        from ..loaders.codon_table import load_codon_table

        return cls(load_codon_table(fp))
