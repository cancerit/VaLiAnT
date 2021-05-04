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
from typing import List, Iterable, Tuple, Optional, FrozenSet
from ..loaders.tsv import load_tsv
from ..models.base import PositionRange
from ..enums import TargetonMutator
from ..utils import parse_mutators


CSV_HEADER = [
    'seq_id',
    'targeton_start',
    'targeton_end',
    'r2_start',
    'r2_end',
    'action_vector'
]


@dataclass(frozen=True)
class CDNATargetonConfig:
    __slots__ = {'seq_id', 'targeton_range', 'r2_range', 'mutators'}

    seq_id: str
    targeton_range: PositionRange
    r2_range: PositionRange
    mutators: FrozenSet[TargetonMutator]

    @classmethod
    def from_row(cls, row: List[str]) -> CDNATargetonConfig:
        return cls(
            row[0],
            PositionRange(int(row[1]), int(row[2])),
            PositionRange(int(row[3]), int(row[4])),
            frozenset(parse_mutators(row[5])))


@dataclass
class CDNATargetonConfigCollection:
    __slots__ = {'cts'}

    cts: List[CDNATargetonConfig]

    def __init__(self, cts: Iterable[CDNATargetonConfig]) -> None:
        self.cts = list(cts)

    def __len__(self) -> int:
        return len(self.cts)

    @property
    def mutators(self) -> FrozenSet[TargetonMutator]:
        return frozenset.union(*[
            ctc.mutators
            for ctc in self.cts
        ])

    @classmethod
    def load(cls, fp: str) -> CDNATargetonConfigCollection:
        return cls(map(CDNATargetonConfig.from_row, load_tsv(fp, CSV_HEADER)))

    @property
    def sequence_ids(self) -> FrozenSet[str]:
        return frozenset(t.seq_id for t in self.cts)
