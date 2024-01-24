########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2023, 2024 Genome Research Ltd
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

from collections import defaultdict
from dataclasses import dataclass
import os
from sqlite3 import Connection
from typing import ClassVar

from pysam import VariantFile

from .annotation import Annotation
from .constants import REVCOMP_OLIGO_NAME_SUFFIX
from .db import cursor
from .enums import MutationType, SrcType
from .experiment_meta import ExperimentMeta
from .loaders.experiment import ExperimentConfig
from .mave_hgvs import get_mave_nt
from .meta_row import MetaRow, sql_select_meta
from .oligo_generation_info import OligoGenerationInfo
from .options import Options
from .seq import Seq
from .sge_config import SGEConfig
from .strings.dna_str import DnaStr
from .strings.strand import Strand
from .utils import bool_to_int_str, reverse_complement
from .variant import Variant
from .vcf_writer import open_vcf_write


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


def _write_field(fh, s: str | None) -> None:
    if s:
        fh.write(s)
    fh.write(',')


def get_transcript_frag(gene_id: str | None, transcript_id: str | None) -> str:
    return (
        f"{transcript_id}.{gene_id}" if transcript_id and gene_id else
        'NO_TRANSCRIPT'
    )


@dataclass(slots=True)
class MetaTable:
    CSV_HEADER: ClassVar[str] = ','.join(META_CSV_FIELDS)

    # TODO: generalise for cDNA support
    cfg: SGEConfig
    exp: ExperimentMeta
    opt: Options
    exp_cfg: ExperimentConfig
    seq: Seq
    alt_seq: Seq
    ppe_mut_types: list[MutationType]
    annot: Annotation | None

    def _write_header(self, fh) -> None:
        fh.write(self.CSV_HEADER)
        fh.write('\n')

    def to_csv(self, conn: Connection, targeton_name: str) -> OligoGenerationInfo:
        def get_path(suffix: str, ext: str) -> str:
            return self.cfg.get_output_file_path(f"{targeton_name}_{suffix}.{ext}")

        def get_csv_path(suffix: str) -> str:
            return get_path(suffix, 'csv')

        def get_vcf_path(suffix: str) -> str:
            return get_path(suffix, 'vcf')

        # Generate file names
        unique_fp = get_csv_path('unique')
        meta_fp = get_csv_path('meta')
        meta_excl_fp = get_csv_path('meta_excluded')
        vcf_ref_fp = get_vcf_path('ref')
        vcf_pam_fp = get_vcf_path('pam')
        # TODO: extra VCF for background?

        # Map: oligo sequence -> oligo names
        unique_oligos = defaultdict(list)

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
        pam_mut_annot = ';'.join([x.value for x in self.ppe_mut_types])

        gene_id = None
        transcript_id = None
        if self.annot:
            transcript_id = self.annot.transcript_id
            gene_id = self.annot.gene_id

        is_rc = self.opt.should_rc(Strand(strand))

        # TODO: generalise for cDNA support
        src_type = SrcType.REF.value

        info = OligoGenerationInfo()

        def get_full_oligo(frag: str) -> str:
            """Add the adaptors to the oligonucleotide sequence"""

            return f"{self.cfg.adaptor_5}{frag}{self.cfg.adaptor_3}"

        def get_mave_hgvs(pos: int, ref: str, alt: str | None = None) -> str:
            """Compile the MAVE-HGVS string"""

            return get_mave_nt(pos, ref_range.start, v.type, ref, alt or mr.alt)

        def get_oligo_name(src: str, v: Variant) -> str:
            tr_frag = get_transcript_frag(gene_id, transcript_id)
            var_frag = v.get_oligo_name_frag()
            rc = REVCOMP_OLIGO_NAME_SUFFIX if is_rc else ''
            return f"{tr_frag}_{contig}:{var_frag}_{src}{rc}"

        def vcf_write(
            vcf: VariantFile,
            start: int,
            end: int,
            ref: str,
            alt: str,
            mutator: str,
            oligo_name: str,
            vcf_alias: str | None,
            vcf_var_id: str | None,
            sge_ref: str | None = None
        ) -> None:

            # Populate INFO field
            vcf_info = {
                'SGE_SRC': mutator,
                'SGE_OLIGO': oligo_name
            }
            if vcf_alias:
                vcf_info['SGE_VCF_ALIAS'] = vcf_alias
                vcf_info['SGE_VCF_VAR_ID'] = vcf_var_id or ''
            if sge_ref:
                vcf_info['SGE_REF'] = sge_ref

            # Write record to file
            vcf.write(vcf.new_record(
                alleles=(ref, alt),
                contig=contig,
                start=start,
                stop=end,
                info=vcf_info))

        with (
            open(meta_fp, 'w') as meta_fh,
            open(meta_excl_fp, 'w') as meta_excl_fh,
            open_vcf_write(vcf_ref_fp, contig) as vcf_ref_fh,
            open_vcf_write(vcf_pam_fp, contig) as vcf_pam_fh
        ):
            self._write_header(meta_fh)
            self._write_header(meta_excl_fh)

            def wf(s: str | None) -> None:
                _write_field(fh, s)

            with cursor(conn) as cur:
                n = cur.execute("select count(*) from v_meta").fetchone()[0]
                print(f"count: {n}")
                it = cur.execute(sql_select_meta)
                while r := it.fetchone():
                    mr = MetaRow(*r)
                    v = mr.to_variant()

                    oligo_no_adapt = mr.oligo
                    if is_rc:
                        oligo_no_adapt = reverse_complement(oligo_no_adapt)
                    oligo = get_full_oligo(oligo_no_adapt)

                    # TODO: correct position for background offsetting
                    ref_pos = v.pos

                    if mr.ref:
                        ref_ref = self.seq.substr(v.ref_range, rel=False)
                        pam_ref = self.alt_seq.substr(mr.ref_range, rel=False)
                    else:
                        ref_ref = ''
                        pam_ref = ''

                    mave_nt_ref = get_mave_hgvs(ref_pos, ref_ref)
                    mave_nt = get_mave_hgvs(v.pos, pam_ref)

                    src = mr.vcf_alias or mr.mutator
                    pam_var = Variant(v.pos, DnaStr(pam_ref), DnaStr(mr.alt))
                    oligo_name = get_oligo_name(src, pam_var)

                    if mr.mutator in {'aa', 'ala', 'stop'}:
                        mr.mutation_type = ''

                    pam_alt = mr.alt

                    if mr.overlaps_codon:
                        pam_range = mr.pam_ref_range
                        if pam_range.start != v.pos or pam_range.end != v.ref_end:

                            # Correct REF and start position in PAM protected codons
                            pam_codon_ref = self.alt_seq.substr(pam_range, rel=False)
                            # Assumption: the targeton start can't change due to background variants
                            # TODO: this may be violated whenever the full transcript is altered,
                            #  unless it is corrected earlier on (once per targeton)
                            oligo_seq = Seq(ref_range.start, DnaStr(mr.oligo))
                            pam_alt = oligo_seq.substr(
                                pam_range.offset_end(v.alt_ref_delta), rel=False)

                            assert len(pam_codon_ref) <= 3
                            mave_nt = get_mave_hgvs(pam_range.start, pam_codon_ref, alt=pam_alt)

                    # Evaluate oligonucleotide length

                    oligo_length = len(oligo)

                    if oligo_length < self.opt.oligo_min_length:
                        info.too_short += 1
                        fh = meta_excl_fh

                    elif oligo_length > self.opt.oligo_max_length:
                        info.too_long += 1
                        fh = meta_excl_fh

                    else:
                        info.in_range += 1
                        fh = meta_fh

                        # Populate unique collection
                        unique_oligos[oligo].append(oligo_name)

                        vcf_start = mr.pos - 1
                        vcf_end = mr.end - 1

                        # TODO: correct insertion and deletion:
                        #  position, REF and ALT prefix

                        # Write REF VCF record
                        vcf_write(
                            vcf_ref_fh,
                            vcf_start,
                            vcf_end,
                            ref_ref or 'N',
                            mr.alt or 'N',
                            mr.mutator,
                            oligo_name,
                            mr.vcf_alias,
                            mr.vcf_var_id)

                        # Write PAM VCF record
                        # TODO: correct position?
                        vcf_write(
                            vcf_pam_fh,
                            vcf_start,
                            vcf_end,
                            pam_ref or 'N',
                            pam_alt or 'N',
                            mr.mutator,
                            oligo_name,
                            mr.vcf_alias,
                            mr.vcf_var_id,
                            sge_ref=ref_ref)

                    # Write metadata table fields
                    #  (either to the default or to the excluded file)

                    # 1. oligo_name
                    wf(oligo_name)

                    # 2. species
                    wf(species)

                    # 3. assembly
                    wf(assembly)

                    # 4. gene_id
                    wf(gene_id)

                    # 5. transcript_id
                    wf(transcript_id)

                    # 6. src_type
                    wf(src_type)

                    # 7. ref_chr
                    wf(contig)

                    # 8. ref_strand
                    wf(strand)

                    # 9. ref_start
                    wf(ref_start)

                    # 10. ref_end
                    wf(ref_end)

                    # 11. revc
                    wf(revc)

                    # 12. ref_seq
                    wf(ref_seq)

                    # 13. pam_seq
                    wf(pam_seq)

                    # 14. vcf_alias
                    wf(mr.vcf_alias)

                    # 15. vcf_var_id
                    wf(mr.vcf_var_id)

                    # 16. mut_position
                    wf(str(mr.pos))

                    # 17. ref
                    wf(pam_ref)

                    # 18. new
                    wf(v.alt)

                    # 19. ref_aa
                    wf(mr.ref_aa)

                    # 20. alt_aa
                    wf(mr.alt_aa)

                    # 21. mut_type
                    wf(mr.mutation_type)

                    # 22. mutator
                    wf(mr.mutator)

                    # 23. oligo_length
                    wf(str(oligo_length))

                    # 24. mseq
                    wf(oligo)

                    # 25. mseq_no_adapt
                    wf(oligo_no_adapt)

                    # 26. pam_mut_annot
                    wf(pam_mut_annot)

                    # 27. pam_mut_sgrna_id
                    wf(mr.sgrna_ids)

                    # 28. mave_nt
                    wf(mave_nt)

                    # 29. mave_nt_ref
                    wf(mave_nt_ref)

                    # 30. vcf_var_in_const
                    fh.write(str(mr.in_const))

                    fh.write('\n')

        if info.out_of_range_n == 0:
            os.unlink(meta_excl_fp)

        if unique_oligos:

            # Write unique sequences
            with open(unique_fp, 'w') as fh:
                fh.write('oligo_name,mseq\n')
                for seq, names in unique_oligos.items():
                    names.sort()
                    fh.write(names[0])
                    fh.write(',')
                    fh.write(seq)
                    fh.write('\n')

        return info
