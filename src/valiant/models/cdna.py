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

from dataclasses import dataclass
from typing import Optional, Sized, Tuple
from .base import PositionRange, StrandedPositionRange
from .codon_table import CodonTable, START_CODON
from .sequences import Sequence
from .sequence_info import SequenceInfo


@dataclass(frozen=True)
class CDNA(Sized):
    __slots__ = {'seq', 'seq_info'}

    seq: Sequence
    seq_info: SequenceInfo

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def range(self) -> StrandedPositionRange:
        return StrandedPositionRange(1, len(self), '+')

    def get_subsequence(self, pr: PositionRange) -> Sequence:
        return self.seq.get_rel_subsequence(pr)

    def get_subsequence_string(self, pr: PositionRange) -> str:
        return self.get_subsequence(pr).sequence

    def get_subsequence_string_before(self, pos: int, length: Optional[int] = None) -> str:
        r = self.range.get_subrange_before(pos, length=length)
        return self.get_subsequence_string(r) if r is not None else ''

    def get_subsequence_string_after(self, pos: int, length: Optional[int] = None) -> str:
        r = self.range.get_subrange_after(pos, length=length)
        return self.get_subsequence_string(r) if r is not None else ''


@dataclass(frozen=True)
class AnnotatedCDNA(CDNA):
    __slots__ = {'seq', 'seq_info', 'cds_range'}

    cds_range: StrandedPositionRange

    def _is_cds_overlap_valid(self, pr: PositionRange) -> bool:
        cds_start: int = self.cds_range.start
        if pr.start < cds_start:
            if pr.end >= cds_start:
                return False
        elif pr.end > self.cds_range.end:
            return False
        return True

    def get_extended_subsequence(self, pr: PositionRange) -> Tuple[Sequence, str, str]:
        if not self._is_cds_overlap_valid(pr):
            raise ValueError("Position range spans CDS and non-CDS regions!")

        # Get subsequence
        seq: Sequence = self.get_subsequence(pr)

        # Get 5' CDS extension
        delta_start: int = pr.start - self.cds_range.start
        ext_5p_len: int = delta_start % 3
        ext_5p: str = self.get_subsequence_string_before(
            pr.start, ext_5p_len) if ext_5p_len != 0 else ''

        # Get 3' CDS extension
        ext_3p_len: int = (3 - (len(pr) + ext_5p_len) % 3) % 3
        ext_3p: str = self.get_subsequence_string_after(
            pr.end, ext_3p_len) if ext_3p_len != 0 else ''

        return seq, ext_5p, ext_3p

    def get_triplet_at(self, pos: int) -> str:
        return self.get_subsequence(PositionRange(pos, pos + 2)).sequence

    def __post_init__(self) -> None:
        n: int = len(self)
        if len(self.cds_range) % 3:
            raise ValueError("Invalid CDS length!")
        if self.cds_range.start > n or self.cds_range.end > n:
            raise ValueError("CDS out of range!")

    # TODO: support non-ATG start codons?
    def validate(self, codon_table: CodonTable) -> None:
        start_codon: str = self.get_triplet_at(self.cds_range.start)
        if start_codon != START_CODON:
            raise ValueError(f"Invalid start codon: {start_codon}!")

        stop_codon: str = self.get_triplet_at(self.cds_range.end - 2)
        if stop_codon not in codon_table.stop_codons:
            raise ValueError(f"Invalid stop codon: {stop_codon}!")
