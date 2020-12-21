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
