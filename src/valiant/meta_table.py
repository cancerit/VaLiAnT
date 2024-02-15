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

import os
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from sqlite3 import Connection
from typing import ClassVar

from .config import BaseConfig
from .constants import REVCOMP_OLIGO_NAME_SUFFIX
from .db import cursor
from .enums import MutationType, SrcType, VariantType
from .experiment_meta import ExperimentMeta
from .genomic_position_offsets import GenomicPositionOffsets
from .mave_hgvs import get_mave_nt
from .meta_row import MetaRow, sql_select_meta
from .oligo_generation_info import OligoGenerationInfo
from .options import Options
from .seq import Seq
from .strings.dna_str import DnaStr
from .strings.strand import Strand
from .transcript_info import TranscriptInfo
from .utils import bool_to_int_str, reverse_complement
from .variant import Variant
from .vcf_writer import Allele, VcfVariant, open_vcf_write, write_vcf_record


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


def get_sge_oligo_name(gene_id: str | None, transcript_id: str | None, contig: str, is_rc: bool, src: str, v: Variant) -> str:
    tr_frag = get_transcript_frag(gene_id, transcript_id)
    var_frag = v.get_oligo_name_frag()
    rc = REVCOMP_OLIGO_NAME_SUFFIX if is_rc else ''
    return f"{tr_frag}_{contig}:{var_frag}_{src}{rc}"


def get_cdna_oligo_name(gene_id: str | None, transcript_id: str | None, seq_id: str, src: str, v: Variant):
    # BRCA1_NP_009225.1_sense_ENST00000357654.9.ENSG00000012048.23_1705_1707_AAA>AAC_aa
    tr_frag = get_transcript_frag(gene_id, transcript_id)
    var_frag = v.get_oligo_name_frag()
    return f"{seq_id}_{tr_frag}_{var_frag}_{src}"


@dataclass(slots=True)
class MetaTable:
    CSV_HEADER: ClassVar[str] = ','.join(META_CSV_FIELDS)

    src_type: SrcType
    seq_id: str | None
    cfg: BaseConfig
    exp: ExperimentMeta
    opt: Options
    gpo: GenomicPositionOffsets | None
    # TODO: support cDNA
    contig: str
    strand: Strand
    seq: Seq
    alt_seq: Seq
    ppe_mut_types: list[MutationType]
    transcript: TranscriptInfo | None

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
        contig = self.contig
        ref_seq = self.seq.s
        ref_range = self.seq.get_range()
        ref_start = str(ref_range.start)
        ref_end = str(ref_range.end)

        gene_id = None
        transcript_id = None
        if self.transcript:
            transcript_id = self.transcript.transcript_id
            gene_id = self.transcript.gene_id

        src_type = self.src_type.value

        info = OligoGenerationInfo()

        assert self.seq.prev_nt or self.seq.start == 1
        # TODO: correct if the altered sequence has global background variants applied upstream
        assert self.alt_seq.prev_nt or self.alt_seq.start == 1

        def get_full_oligo(frag: str) -> str:
            """Add the adaptors to the oligonucleotide sequence"""

            return f"{self.cfg.adaptor_5}{frag}{self.cfg.adaptor_3}"

        def get_mave_hgvs(pos: int, t: VariantType, ref: str, alt: str | None = None) -> str:
            """Compile the MAVE-HGVS string"""

            return get_mave_nt(pos, ref_range.start, t, ref, alt or mr.alt)

        match self.src_type:
            case SrcType.REF:
                strand = self.strand
                is_rc = self.opt.should_rc(Strand(strand))
                get_oligo_name_f = partial(
                    get_sge_oligo_name, gene_id, transcript_id, contig, is_rc)
                # tr_frag = get_transcript_frag(gene_id, transcript_id)
                # var_frag = v.get_oligo_name_frag()
                # rc = REVCOMP_OLIGO_NAME_SUFFIX if is_rc else ''
                # return f"{tr_frag}_{contig}:{var_frag}_{src}{rc}"
                pam_seq = self.alt_seq.s
                pam_mut_annot = ';'.join([x.value for x in self.ppe_mut_types])

            case SrcType.CDNA:
                assert self.seq_id
                strand = ''
                is_rc = False
                get_oligo_name_f = partial(get_cdna_oligo_name, gene_id, transcript_id, self.seq_id)
                pam_seq = None
                pam_mut_annot = None

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
                it = cur.execute(sql_select_meta)
                while r := it.fetchone():
                    mr = MetaRow(*r)
                    ref_var = mr.to_ref_variant()
                    alt_var = mr.to_alt_variant()
                    var_type = ref_var.type

                    oligo_no_adapt = mr.oligo
                    if is_rc:
                        oligo_no_adapt = reverse_complement(oligo_no_adapt)
                    oligo = get_full_oligo(oligo_no_adapt)

                    if mr.ref:
                        ref_ref = self.seq.substr(ref_var.ref_range, rel=False)
                        pam_ref = self.alt_seq.substr(mr.alt_ref_range, rel=False)
                    else:
                        ref_ref = ''
                        pam_ref = ''

                    mave_pos = ref_var.pos
                    mave_nt_ref = get_mave_hgvs(mave_pos, var_type, ref_ref)
                    mave_nt = get_mave_hgvs(mave_pos, var_type, pam_ref)  # ALT position instead?

                    src = mr.vcf_alias or mr.mutator
                    pam_var = Variant(ref_var.pos, DnaStr(pam_ref), DnaStr(mr.alt))
                    oligo_name = get_oligo_name_f(src, pam_var)

                    # TODO: consider dropping (just for backwards compatibility)
                    if mr.mutator in {'aa', 'ala', 'stop'}:
                        mr.mutation_type = ''

                    pam_ref_start = ref_var.pos

                    # Default SGE_REF (assuming no PPE interaction)
                    ref_ref_all = Allele(ref_ref)
                    sge_ref_all = ref_ref_all.clone()
                    pam_ref_all = Allele(pam_ref)
                    ref_alt_all = Allele(mr.alt)
                    pam_alt_all = ref_alt_all.clone()

                    add_vcf_nt = var_type != VariantType.SUBSTITUTION

                    # This will only be set for custom variants, by definition being from reference
                    ref_nt: str | None = mr.vcf_nt or None
                    alt_nt: str | None = ref_nt

                    if add_vcf_nt:
                        ref_ref_all.nt = ref_nt
                        sge_ref_all.nt = ref_nt
                        pam_ref_all.nt = ref_nt
                        ref_alt_all.nt = ref_nt
                        pam_alt_all.nt = ref_nt

                    if (
                        add_vcf_nt and
                        not ref_nt and
                        mr.mutator != 'custom'
                    ):
                        ref_nt = self.seq.get_nt(ref_var.pos - 1)
                        alt_nt = self.alt_seq.get_nt(alt_var.pos - 1)
                        sge_ref_all.nt = ref_nt
                        ref_ref_all.nt = ref_nt
                        pam_ref_all.nt = alt_nt
                        ref_alt_all.nt = alt_nt
                        pam_alt_all.nt = alt_nt

                    if mr.overlaps_codon:
                        # This range may exceed the span of a codon if so does the mutation
                        pam_range = mr.pam_ref_range
                        if pam_range.start != alt_var.pos or pam_range.end != alt_var.ref_end:

                            # Correct REF and start position in PAM protected codons
                            pam_codon_ref = self.alt_seq.substr(pam_range, rel=False)
                            codon_ref = self.seq.substr(pam_range, rel=False)
                            # Assumption: the targeton start can't change due to background variants
                            # TODO: this may be violated whenever the full transcript is altered,
                            #  unless it is corrected earlier on (once per targeton)
                            oligo_seq = Seq(ref_range.start, DnaStr(mr.oligo))
                            pam_alt = oligo_seq.substr(
                                pam_range.offset_end(alt_var.alt_ref_delta), rel=False)

                            sge_ref_all.s = codon_ref
                            pam_ref_all.s = pam_codon_ref
                            pam_alt_all.s = pam_alt

                            if pam_range.start <= alt_var.pos - 1:
                                sge_ref_all.nt = None
                                pam_ref_all.nt = None
                                pam_alt_all.nt = None

                            # Lift over to reference
                            pam_ref_start = (
                                self.gpo.alt_to_ref_position(pam_range.start) if self.gpo else
                                pam_range.start
                            )
                            assert pam_ref_start is not None

                            mave_nt = get_mave_hgvs(
                                pam_ref_start, var_type, pam_codon_ref, alt=pam_alt)

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

                        if self.src_type == SrcType.REF:

                            # TODO: verify it is fine to drop the stored end position
                            vcf_ref = VcfVariant.from_partial(
                                ref_var.pos - 1,
                                ref_ref_all,
                                ref_alt_all)

                            vcf_pam = VcfVariant.from_partial(
                                pam_ref_start - 1,
                                pam_ref_all,
                                pam_alt_all,
                                sge_ref=sge_ref_all)

                            # Write REF VCF record
                            write_vcf_record(
                                vcf_ref_fh,
                                contig,
                                vcf_ref,
                                mr.mutator,
                                oligo_name,
                                mr.vcf_alias,
                                mr.vcf_var_id)

                            # Write PAM VCF record
                            write_vcf_record(
                                vcf_pam_fh,
                                contig,
                                vcf_pam,
                                mr.mutator,
                                oligo_name,
                                mr.vcf_alias,
                                mr.vcf_var_id)

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
                    wf(str(ref_var.pos))

                    # 17. ref
                    wf(pam_ref)

                    # 18. new
                    wf(ref_var.alt)

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

        if info.in_range == 0:
            os.unlink(meta_fp)

        if info.in_range == 0 or self.src_type == SrcType.CDNA:
            os.unlink(vcf_ref_fp)
            os.unlink(vcf_pam_fp)

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
