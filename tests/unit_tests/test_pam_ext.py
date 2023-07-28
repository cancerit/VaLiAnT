########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023 Genome Research Ltd
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

import numpy as np
import pandas as pd
import pytest

from valiant.constants import META_PAM_MUT_SGRNA_ID, META_REF, META_NEW, META_MUT_POSITION, META_REF_START, META_REF_SEQ, META_PAM_SEQ, META_MSEQ_NO_ADAPT_NO_RC, META_PAM_CODON_REF, META_PAM_CODON_ALT
from valiant.metadata_utils import set_pam_extended_ref_alt, set_ref_meta
from valiant.models.base import GenomicRange, GenomicPosition
from valiant.models.new_pam import CdsPamBgAltSeqBuilder
from valiant.models.pam_protection import PamVariant
from valiant.models.targeton import PamProtCDSTargeton
from valiant.models.variant import DeletionVariant


def _get_test_pos(pos: int):
    return GenomicPosition('X', pos)


def _get_test_range(start: int, length: int):
    return GenomicRange('X', start, start + length - 1, '+')


def seq_to_meta(seq: CdsPamBgAltSeqBuilder, mut: DeletionVariant):
    return {
        META_REF_START: seq.start,
        META_REF_SEQ: seq.ref_seq,
        META_PAM_SEQ: seq.alt_seq,
        META_MUT_POSITION: mut.genomic_position.position,
        META_REF: mut.ref,
        META_NEW: None,
        META_MSEQ_NO_ADAPT_NO_RC: mut.mutate(seq.alt_seq, seq.start, ref_check=False),
        META_PAM_MUT_SGRNA_ID: ','.join(x.sgrna_id for x in seq.pam_variants)  # format?
    }


# Mutated first protected PAM codon
# TAC[AAA]GAC
#   ^       ^
#  PAM     PAM
@pytest.mark.parametrize('del_start,del_length,pam_ref_exp,pam_alt_exp', [
    (100, 1, 'TAC', 'AC'),
    (101, 1, 'AC', 'C'),
    (102, 1, 'C', ''),
    (127, 1, 'GAC', 'AC'),
    (127, 2, 'GAC', 'C'),
    (128, 1, 'AC', 'C'),
    (129, 1, 'C', ''),
    (101, 27, 'ACAAAAAAAAAAAAAAAAAAAAAAAAGAC', 'AC'),
    (101, 28, 'ACAAAAAAAAAAAAAAAAAAAAAAAAGAC', 'C'),
    (120, 8, 'AAAAAAAGAC', 'AC')
])
def test_set_pam_extended_ref_alt(del_start, del_length, pam_ref_exp, pam_alt_exp):
    assert del_length > 0

    df = pd.DataFrame()
    sgrna_id = 'sgrna-a'
    start = 100
    assert del_start >= start

    seq_length = 30
    ref = 'T' + 'A' * (seq_length - 4) + 'GAA'
    gr = _get_test_range(start, seq_length)

    # Set PAM protection variants
    pam_variant_a = PamVariant(_get_test_pos(gr.start + 2), 'A', 'C', sgrna_id)
    pam_variant_b = PamVariant(_get_test_pos(gr.end), 'A', 'C', sgrna_id)

    # Apply PAM protection variants to the reference sequence
    seq = CdsPamBgAltSeqBuilder.from_ref(gr, ref, [], [
        pam_variant_a,
        pam_variant_b
    ])
    targeton = PamProtCDSTargeton(seq)

    mut = DeletionVariant(_get_test_pos(del_start), 'G' * del_length)

    # Add mutation to the metadata table
    df = df.append(seq_to_meta(seq, mut), ignore_index=True)

    # Cast integer fields
    for f in [
        META_REF_START,
        META_MUT_POSITION
    ]:
        df[f] = df[f].astype(np.int64)

    # Add metadata
    set_ref_meta(df)
    set_pam_extended_ref_alt(df, [targeton])

    # Check reference and alternative including the protected PAM codon
    assert df.loc[0, META_PAM_CODON_REF] == pam_ref_exp
    assert df.loc[0, META_PAM_CODON_ALT] == pam_alt_exp
