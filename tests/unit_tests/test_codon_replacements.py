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

import numpy as np
import pytest
from valiant.enums import VariantType
from valiant.models.base import GenomicRange
from valiant.models.new_pam import CdsPamBgAltSeqBuilder
from valiant.models.pam_protected_reference_sequence import PamProtectedReferenceSequence
from valiant.models.sequences import ReferenceSequence
from valiant.models.targeton import CDSTargeton
from .utils import get_aux_tables, seq2triplets

CONST_CODON_METHODS = {
    'A': 'get_ala_mutations',
    'STOP': 'get_stop_mutations'
}

aux = get_aux_tables()
sub_var_type = np.int8(VariantType.SUBSTITUTION.value)


@pytest.mark.parametrize('aa,strand,seq,exp_mseq', [
    ('A', '+', 'GCAAAATTT', ['GCCAAATTT', 'GCAGCCTTT', 'GCAAAAGCC']),
    ('A', '+', 'GCCAAATTT', ['GCCGCCTTT', 'GCCAAAGCC']),
    ('A', '-', 'TGCAAATTT', ['GGCAAATTT', 'TGCGGCTTT', 'TGCAAAGGC']),
    ('A', '+', 'GATTTT', ['GCCTTT', 'GATGCC']),
    ('STOP', '+', 'TAACCCGGG', ['TGACCCGGG', 'TAATGAGGG', 'TAACCCTGA']),
    ('STOP', '-', 'TTACCCGGG', ['TCACCCGGG', 'TTATCAGGG', 'TTACCCTCA'])
])
def test_get_snvre_aa_mutations(aa, strand, seq, exp_mseq):

    # Generate target
    gr = GenomicRange('X', 10, 10 + len(seq) - 1, strand)
    t = CDSTargeton(CdsPamBgAltSeqBuilder.from_ref(gr, seq, [], []))

    # Generate codon substitution mutations
    mc = getattr(t, CONST_CODON_METHODS[aa])(aux_tables=aux)

    # Check results
    assert mc.df.alt_aa.unique()[0] == aa
    assert np.array_equal(mc.df.mseq.to_numpy(), np.array(exp_mseq))
    assert mc.df.var_type.unique()[0] == sub_var_type
    assert np.array_equal(
        mc.df.ref.astype('string').to_numpy(),
        np.array([
            triplet
            for i, triplet in enumerate(seq2triplets(seq))
            if 3 * i in mc.df.mut_position.unique()
        ]))
