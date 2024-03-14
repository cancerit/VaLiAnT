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

from contextlib import contextmanager, nullcontext
from io import TextIOWrapper
import os
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from sqlite3 import Connection
from typing import ClassVar, Generator

from pysam import VariantFile

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
    'vcf_var_in_const',
    'background_variants',
    'background_seq'
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


def get_sge_oligo_no_op_name(gene_id: str | None, transcript_id: str | None, contig: str, is_rc: bool) -> str:
    tr_frag = get_transcript_frag(gene_id, transcript_id)
    rc = REVCOMP_OLIGO_NAME_SUFFIX if is_rc else ''
    return f"{tr_frag}_{contig}_no_op{rc}"


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


def write_no_op_meta_record(
    fh: TextIOWrapper,
    oligo_name: str,
    species: str,
    assembly: str,
    gene_id: str | None,
    transcript_id: str | None,
    contig: str,
    strand: str,
    ref_start: str,
    ref_end: str,
    revc: str,
    ref_seq: str,
    pam_seq: str | None,
    oligo_length: int,
    oligo: str,
    oligo_no_adapt: str,
    pam_mut_annot: str | None,
    background_variants: str,
    background_seq: str
) -> None:
    write_meta_record(
        fh,
        oligo_name,
        species,
        assembly,
        gene_id,
        transcript_id,
        '',
        contig,
        strand,
        ref_start,
        ref_end,
        revc,
        ref_seq,
        pam_seq,
        '',
        '',
        -1,
        '',
        '',
        '',
        '',
        '',
        '',
        oligo_length,
        oligo,
        oligo_no_adapt,
        pam_mut_annot,
        '',
        '',
        '',
        0,
        background_variants,
        background_seq)


def write_meta_record(
    fh: TextIOWrapper,
    oligo_name: str,
    species: str,
    assembly: str,
    gene_id: str | None,
    transcript_id: str | None,
    src_type: str,
    contig: str,
    strand: str,
    ref_start: str,
    ref_end: str,
    revc: str,
    ref_seq: str,
    pam_seq: str | None,
    vcf_alias: str,
    vcf_var_id: str,
    ref_var_pos: int,
    pam_ref: str,
    ref_var_alt: str,
    ref_aa: str,
    alt_aa: str,
    mutation_type: str,
    mutator: str,
    oligo_length: int,
    oligo: str,
    oligo_no_adapt: str,
    pam_mut_annot: str | None,
    sgrna_ids: str,
    mave_nt: str,
    mave_nt_ref: str,
    in_const: int,
    background_variants: str,
    background_seq: str
) -> None:
    # 1. oligo_name
    _write_field(fh, oligo_name)

    # 2. species
    _write_field(fh, species)

    # 3. assembly
    _write_field(fh, assembly)

    # 4. gene_id
    _write_field(fh, gene_id)

    # 5. transcript_id
    _write_field(fh, transcript_id)

    # 6. src_type
    _write_field(fh, src_type)

    # 7. ref_chr
    _write_field(fh, contig)

    # 8. ref_strand
    _write_field(fh, strand)

    # 9. ref_start
    _write_field(fh, ref_start)

    # 10. ref_end
    _write_field(fh, ref_end)

    # 11. revc
    _write_field(fh, revc)

    # 12. ref_seq
    _write_field(fh, ref_seq)

    # 13. pam_seq
    _write_field(fh, pam_seq)

    # 14. vcf_alias
    _write_field(fh, vcf_alias)

    # 15. vcf_var_id
    _write_field(fh, vcf_var_id)

    # 16. mut_position
    _write_field(fh, str(ref_var_pos))

    # 17. ref
    _write_field(fh, pam_ref)

    # 18. new
    _write_field(fh, ref_var_alt)

    # 19. ref_aa
    _write_field(fh, ref_aa)

    # 20. alt_aa
    _write_field(fh, alt_aa)

    # 21. mut_type
    _write_field(fh, mutation_type)

    # 22. mutator
    _write_field(fh, mutator)

    # 23. oligo_length
    _write_field(fh, str(oligo_length))

    # 24. mseq
    _write_field(fh, oligo)

    # 25. mseq_no_adapt
    _write_field(fh, oligo_no_adapt)

    # 26. pam_mut_annot
    _write_field(fh, pam_mut_annot)

    # 27. pam_mut_sgrna_id
    _write_field(fh, sgrna_ids)

    # 28. mave_nt
    _write_field(fh, mave_nt)

    # 29. mave_nt_ref
    _write_field(fh, mave_nt_ref)

    # 30. vcf_var_in_const
    _write_field(fh, str(in_const))

    # 31. background_variants
    _write_field(fh, background_variants)

    # 32. background_seq
    fh.write(background_seq)
    fh.write('\n')


@dataclass(slots=True)
class MetaTable:
    CSV_HEADER: ClassVar[str] = ','.join(META_CSV_FIELDS)

    src_type: SrcType
    seq_id: str | None
    cfg: BaseConfig
    exp: ExperimentMeta
    opt: Options
    gpo: GenomicPositionOffsets | None
    contig: str
    strand: Strand
    seq: Seq
    bg_seq: Seq
    alt_seq: Seq
    bg_variants: list[Variant]
    ppe_mut_types: list[MutationType]
    transcript: TranscriptInfo | None

    @contextmanager
    def open_vcf_write(self, fp: str) -> Generator[VariantFile | None, None, None]:
        ctx = (
            open_vcf_write(fp, self.contig) if self.src_type == SrcType.REF else
            nullcontext()
        )
        with ctx as fh:
            yield fh

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

        adaptor_5 = self.cfg.adaptor_5 or ''
        adaptor_3 = self.cfg.adaptor_3 or ''

        gene_id = None
        transcript_id = None
        if self.transcript:
            transcript_id = self.transcript.transcript_id
            gene_id = self.transcript.gene_id

        src_type = self.src_type.value

        # Targeton-level background-specific metadata
        if self.bg_variants:
            background_seq = self.bg_seq.s
            # Coordinates expected to be in the reference coordinate system
            background_variants = ';'.join(
                get_mave_nt(
                    x.pos,
                    ref_range.start,
                    x.type,
                    x.ref or None,
                    x.alt or None
                ) for x in self.bg_variants
            )
        else:
            background_seq = ref_seq
            background_variants = ''

        info = OligoGenerationInfo()

        assert self.seq.prev_nt or self.seq.start == 1
        # TODO: correct if the altered sequence has global background variants applied upstream
        assert self.alt_seq.prev_nt or self.alt_seq.start == 1

        def get_full_oligo(frag: str) -> str:
            """Add the adaptors to the oligonucleotide sequence"""

            return f"{adaptor_5}{frag}{adaptor_3}"

        def get_mave_hgvs(pos: int, t: VariantType, ref: str, alt: str | None = None) -> str:
            """Compile the MAVE-HGVS string"""

            return get_mave_nt(pos, ref_range.start, t, ref, alt or mr.alt)

        match self.src_type:
            case SrcType.REF:
                strand = self.strand
                is_rc = self.opt.should_rc(Strand(strand))
                get_oligo_name_f = partial(
                    get_sge_oligo_name, gene_id, transcript_id, contig, is_rc)
                pam_seq = self.alt_seq.s
                pam_mut_annot = ';'.join([x.value for x in self.ppe_mut_types])

            case SrcType.CDNA:
                assert self.seq_id
                strand = ''
                is_rc = False
                get_oligo_name_f = partial(
                    get_cdna_oligo_name, gene_id, transcript_id, self.seq_id)
                pam_seq = None
                pam_mut_annot = None

        # Name of the no-op oligonucleotide
        no_op_oligo_name = get_sge_oligo_no_op_name(gene_id, transcript_id, contig, is_rc)

        with (
            open(meta_fp, 'w') as meta_fh,
            open(meta_excl_fp, 'w') as meta_excl_fh,
            self.open_vcf_write(vcf_ref_fp) as vcf_ref_fh,
            self.open_vcf_write(vcf_pam_fp) as vcf_pam_fh
        ):
            self._write_header(meta_fh)
            self._write_header(meta_excl_fh)

            with cursor(conn) as cur:

                if self.opt.include_no_op_oligo:

                    # Prepare no-op oligonucleotide sequences
                    oligo_no_adapt = background_seq
                    if is_rc:
                        oligo_no_adapt = reverse_complement(oligo_no_adapt)
                    oligo = get_full_oligo(oligo_no_adapt)
                    oligo_length = len(oligo)

                    # Evalute oligonucleotide length
                    if info.eval_in_range(self.opt, oligo_length):
                        fh = meta_fh
                        unique_oligos[oligo].append(no_op_oligo_name)
                    else:
                        fh = meta_excl_fh

                    write_no_op_meta_record(
                        fh,
                        no_op_oligo_name,
                        species,
                        assembly,
                        gene_id,
                        transcript_id,
                        contig,
                        strand,
                        ref_start,
                        ref_end,
                        revc,
                        ref_seq,
                        pam_seq,
                        oligo_length,
                        oligo,
                        oligo_no_adapt,
                        pam_mut_annot,
                        background_variants,
                        background_seq)

                it = cur.execute(sql_select_meta)
                while r := it.fetchone():
                    mr = MetaRow(*r)
                    ref_var = mr.to_ref_variant()
                    alt_var = mr.to_alt_variant()
                    var_type = ref_var.type

                    # Prepare oligonucleotide sequences
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
                    if info.eval_in_range(self.opt, oligo_length):
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
                                vcf_ref_fh,  # type: ignore
                                contig,
                                vcf_ref,
                                mr.mutator,
                                oligo_name,
                                mr.vcf_alias,
                                mr.vcf_var_id)

                            # Write PAM VCF record
                            write_vcf_record(
                                vcf_pam_fh,  # type: ignore
                                contig,
                                vcf_pam,
                                mr.mutator,
                                oligo_name,
                                mr.vcf_alias,
                                mr.vcf_var_id)
                    else:
                        fh = meta_excl_fh

                    # Write metadata table fields
                    #  (either to the default or to the excluded file)

                    write_meta_record(
                        fh,
                        oligo_name,
                        species,
                        assembly,
                        gene_id,
                        transcript_id,
                        src_type, contig,
                        strand,
                        ref_start,
                        ref_end,
                        revc,
                        ref_seq,
                        pam_seq,
                        mr.vcf_alias,
                        mr.vcf_var_id,
                        ref_var.pos,
                        pam_ref,
                        ref_var.alt,
                        mr.ref_aa,
                        mr.alt_aa,
                        mr.mutation_type,
                        mr.mutator,
                        oligo_length,
                        oligo,
                        oligo_no_adapt,
                        pam_mut_annot,
                        mr.sgrna_ids,
                        mave_nt,
                        mave_nt_ref,
                        mr.in_const,
                        background_variants,
                        background_seq)

        if info.out_of_range_n == 0:
            os.unlink(meta_excl_fp)

        if info.in_range == 0:
            os.unlink(meta_fp)

        if info.in_range == 0 and self.src_type == SrcType.REF:
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
