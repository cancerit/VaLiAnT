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

from functools import lru_cache
import os
import pathlib
import re
from typing import List, Type, Tuple
import numpy as np
import pandas as pd

dna_complement_tr_table = str.maketrans('ACGT', 'TGCA')
dna_re = re.compile('^[ACGT]+$')

i8_0: np.int8 = np.int8(0)
i8_1: np.int8 = np.int8(1)
i8_2: np.int8 = np.int8(2)

i8_1_2: List[np.int8] = [i8_1, i8_2]
i8_2_1: List[np.int8] = [i8_2, i8_1]


def get_data_file_path(fp: str) -> str:
    return os.path.join(pathlib.Path(__file__).parent.absolute(), 'data', fp)


def is_dna(seq: str) -> bool:
    return dna_re.match(seq) is not None


def is_strand(strand: str) -> bool:
    return strand == '+' or strand == '-'


def validate_strand(strand: str) -> None:
    if not is_strand(strand):
        raise ValueError(f"Invalid strand '{strand}'!")


def get_region(chromosome: str, start: int, end: int) -> str:
    return f"{chromosome}:{start}-{end}"


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(dna_complement_tr_table)


def parse_list(s: str) -> List[str]:
    return [
        item for item in [
            raw.strip()
            for raw in s.split(',')
        ]
        if item
    ]


def get_smallest_int_type(max_value: int) -> Type:
    if max_value < 0:
        raise ValueError("Only non-negative integers supported!")
    return (
        np.int8 if max_value <= np.iinfo(np.int8).max else
        np.int16 if max_value <= np.iinfo(np.int16).max else
        np.int32 if max_value <= np.iinfo(np.int32).max else
        np.int64
    )


def get_id_column(rown: int) -> np.ndarray:
    return np.arange(rown, dtype=get_smallest_int_type(rown))


def get_constant_category(s: str, n: int, categories: List[str] = None) -> pd.Categorical:
    return pd.Categorical([s], categories=categories or [s]).repeat(n)


@lru_cache(maxsize=4)
def get_out_of_frame_offset(cds_ext_length: int) -> int:
    return (3 - cds_ext_length) if cds_ext_length > 0 else 0


def get_frame_complement(frame: pd.Series) -> pd.Series:
    return frame.replace(i8_1_2, i8_2_1)


def get_inner_cds_relative_boundaries(seq_len: int, frame: int) -> Tuple[int, int]:
    inner_cds_start: int = 0 if frame == 0 else (3 - frame)
    inner_cds_end: int = seq_len - (seq_len - inner_cds_start) % 3
    return inner_cds_start, inner_cds_end


def get_var_types(var_types: pd.Series) -> List[int]:
    return [
        int(k)
        for k, v in dict(var_types.value_counts()).items()
        if v > 0
    ]
