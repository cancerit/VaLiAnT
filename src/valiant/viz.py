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


def _pad(pad: int, s: str) -> str:
    return s + ' ' * (pad - len(s)) if pad and len(s) < pad else s


def _draw_region(start: int, label: str, gr: GenomicRange, pad: int = 0):
    return _pad(pad, label) + _draw_range(start, gr)


def get_text_diagram(b: PamBgAltSeqBuilder, rsr_ref: ReferenceSequenceRanges, rsr_alt: ReferenceSequenceRanges) -> str:
    # TODO: add CDS annotation (with extensions?)
    # TODO: add ALT sequence (with alignment?)

    s = []

    for seq, rsr in [(b.ref_seq, rsr_ref), (b.alt_seq, rsr_alt)]:
        start: int = rsr.ref_range.start
        pad = 4
        s.append((' ' * pad) + seq)

        bg_variants = b.bg_variants
        if bg_variants:
            bg = {z.start - start: z for z in bg_variants}
            ps = _pad(pad, "BG")
            for i in range(len(seq)):
                ps += ('!' if bg[i].is_frame_shifting else '^') if i in bg else ' '
            s.append(ps)

        pam_variants = b.pam_variants
        if pam_variants:
            pam = {z.start - start for z in pam_variants}
            ps = _pad(pad, "PAM")
            for i in range(len(seq)):
                ps += '^' if i in pam else ' '
            s.append(ps)

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
