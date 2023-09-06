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

from .models.base import GenomicRange
from .models.new_pam import PamBgAltSeqBuilder
from .models.refseq_ranges import ReferenceSequenceRanges


def _draw_range(start: int, gr: GenomicRange) -> str:
    return ' ' * (gr.start - start) + '*' * len(gr)


def _draw_region(start: int, label: str, gr: GenomicRange, pad: int = 0):
    if pad and len(label) < pad:
        label += ' ' * (pad - len(label))
    return label + _draw_range(start, gr)


def get_text_diagram(b: PamBgAltSeqBuilder, rsr: ReferenceSequenceRanges) -> str:
    # TODO: add CDS annotation (with extensions?)
    # TODO: add ALT sequence (with alignment?)
    start: int = rsr.ref_range.start
    pad = 4
    s = [(' ' * pad) + b.ref_seq]

    f = lambda x, y: _draw_region(start, x, y, pad=4)

    if rsr.const_region_1:
        s.append(f('C1', rsr.const_region_1))
    if rsr.target_ragion_1:
        s.append(f('R1', rsr.target_ragion_1))
    s.append(f('R2', rsr.target_ragion_2))
    if rsr.target_ragion_3:
        s.append(f('R3', rsr.target_ragion_3))
    if rsr.const_region_2:
        s.append(f('C2', rsr.const_region_2))

    return '\n'.join(s)
