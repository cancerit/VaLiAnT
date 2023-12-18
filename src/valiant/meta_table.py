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

from dataclasses import dataclass
from sqlite3 import Connection
from typing import ClassVar

from .db import cursor
from .experiment_meta import ExperimentMeta
from .loaders.experiment import ExperimentConfig
from .options import Options
from .seq import Seq
from .sge_config import SGEConfig
from .utils import bool_to_int_str


META_CSV_FIELDS = [
    'oligo_name',
    'species',
    'assembly',
    'gene_id',
    'transcript_id',
    'src_type',
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'revc',
    'ref_seq',
    'pam_seq',
    'vcf_alias',
    'vcf_var_id',
    'mut_position',
    'ref',
    'new',
    'ref_aa',
    'alt_aa',
    'mut_type',
    'mutator',
    'oligo_length',
    'mseq',
    'mseq_no_adapt',
    'pam_mut_annot',
    'pam_mut_sgrna_id',
    'mave_nt',
    'mave_nt_ref',
    'vcf_var_in_const'
]


SQL_SELECT = """
select
    ref_start,
    ref,
    alt,
    ref_aa,
    alt_aa,
    vcf_var_id,
    vcf_alias,
    mutator,
    in_const,
    oligo,
    sgrna_ids
from v_meta
"""


def _write_field(fh, s: str | None) -> None:
    if s:
        fh.write(s)
    fh.write(',')


@dataclass(slots=True)
class MetaRow:
    pos: int
    ref: str
    alt: str
    ref_aa: str
    alt_aa: str
    vcf_var_id: str
    vcf_alias: str
    mutator: str
    in_const: int
    oligo: str
    sgrna_ids: str


@dataclass
class MetaTable:
    CSV_HEADER: ClassVar[str] = ','.join(META_CSV_FIELDS)

    # TODO: generalise for cDNA support
    cfg: SGEConfig
    exp: ExperimentMeta
    opt: Options
    exp_cfg: ExperimentConfig
    seq: Seq
    alt_seq: Seq

    def _write_header(self, fh) -> None:
        fh.write(self.CSV_HEADER)
        fh.write('\n')

    def to_csv(self, conn: Connection, fp: str) -> None:
        species = self.cfg.species
        assembly = self.cfg.assembly
        revc = bool_to_int_str(self.opt.revcomp_minus_strand)
        contig = self.exp_cfg.contig
        strand = self.exp_cfg.strand
        ref_seq = self.seq.s
        ref_range = self.seq.get_range()
        ref_start = str(ref_range.start)
        ref_end = str(ref_range.end)
        pam_seq = self.alt_seq.s

        # TODO: generalise for cDNA support
        src_type = 'ref'

        # TODO: recover proper struct for reporting
        too_short = 0
        too_long = 0

        def get_full_oligo(mr: MetaRow) -> str:
            return f"{self.cfg.adaptor_5}{mr.oligo}{self.cfg.adaptor_3}"

        with open(fp, 'w') as fh:
            self._write_header(fh)

            def wf(s: str | None) -> None:
                _write_field(fh, s)

            with cursor(conn) as cur:
                it = cur.execute(SQL_SELECT)
                while r := it.fetchone():
                    mr = MetaRow(*r)

                    oligo = get_full_oligo(mr)
                    oligo_length = len(oligo)

                    if oligo_length < self.opt.oligo_min_length:
                        too_short += 1
                        continue

                    if oligo_length > self.opt.oligo_max_length:
                        too_long += 1
                        continue

                    # Write fields

                    # oligo_name
                    # TODO
                    wf('<NAME>')

                    # species
                    wf(species)

                    # assembly
                    wf(assembly)

                    # gene_id
                    # TODO
                    wf('<GID>')

                    # transcript_id
                    # TODO
                    wf('<TID>')

                    # src_type
                    wf(src_type)

                    # ref_chr
                    wf(contig)

                    # ref_strand
                    wf(strand)

                    # ref_start
                    wf(ref_start)

                    # ref_end
                    wf(ref_end)

                    # revc
                    wf(revc)

                    # ref_seq
                    wf(ref_seq)

                    # pam_seq
                    wf(pam_seq)

                    # vcf_alias
                    wf(mr.vcf_alias)

                    # vcf_var_id
                    wf(mr.vcf_var_id)

                    # mut_position
                    wf(str(mr.pos))

                    # ref
                    wf(mr.ref)

                    # new
                    wf(mr.alt)

                    # ref_aa
                    wf(mr.ref_aa)

                    # alt_aa
                    wf(mr.alt_aa)

                    # mut_type
                    # TODO
                    wf('<MT>')

                    # mutator
                    wf(mr.mutator)

                    # oligo_length
                    wf(str(oligo_length))

                    # mseq
                    wf(oligo)

                    # mseq_no_adapt
                    wf(mr.oligo)

                    # pam_mut_annot
                    # TODO
                    wf('<PMA>')

                    # pam_mut_sgrna_id
                    wf(mr.sgrna_ids)

                    # mave_nt
                    # TODO
                    wf('<MV1>')

                    # mave_nt_ref
                    # TODO
                    wf('<MV2>')

                    # vcf_var_in_const
                    fh.write(str(mr.in_const))

                    fh.write('\n')
