########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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

from hashlib import md5
import os
import pathlib
import pandas as pd
from valiant.models.base import GenomicRange
from valiant.models.codon_table import CodonTable
from valiant.models.pam_protection import PamProtectedReferenceSequence
from valiant.models.sequences import ReferenceSequence
from valiant.models.snv_table import AuxiliaryTables
from valiant.models.targeton import Targeton
from .constants import DUMMY_PAM_PROTECTION_NT, CODON_TABLE_FP, FRAMES, STRANDS

cds_extension_slices = {
    (a, b): slice(a, -b if b > 0 else None)
    for a in FRAMES
    for b in FRAMES
}


def get_data_file_path(fp):
    return os.path.join(pathlib.Path(__file__).parent.absolute(), 'data', fp)


def get_dummy_pam_protected(seq):
    return DUMMY_PAM_PROTECTION_NT * len(seq)


def get_pam_protected_sequence(seq, pam_protection, chromosome='X', strand='+', pos=1):
    gr = GenomicRange(chromosome, pos, len(seq), strand)
    ref_seq = ReferenceSequence(seq, gr)
    return PamProtectedReferenceSequence.from_reference_sequence(
        ref_seq, get_dummy_pam_protected(seq) if pam_protection else seq)


def get_targeton(seq, pam_protection, chromosome='X', strand='+', pos=1):
    pam_ref_seq = get_pam_protected_sequence(
        seq, pam_protection, chromosome=chromosome, strand=strand, pos=pos)
    return Targeton.from_pam_seq(pam_ref_seq)


def load_codon_table():
    return CodonTable.load(get_data_file_path(CODON_TABLE_FP))


def get_aux_tables(strands=None):
    strands_ = strands if strands is not None else STRANDS
    return AuxiliaryTables(load_codon_table(), strands_, True, True, True)


def load_csv(fp):
    return pd.read_csv(get_data_file_path(fp))


def md5sum(s):
    return md5(s.encode('utf-8')).hexdigest()


def trim_cds_extensions(cds_seq, cds_prefix_length, cds_suffix_length):
    assert len(cds_seq) % 3 == 0
    return cds_seq[cds_extension_slices[cds_prefix_length, cds_suffix_length]]


def seq2triplets(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]
