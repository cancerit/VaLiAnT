########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023, 2024 Genome Research Ltd
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

import logging
import sys
from functools import partial
from sqlite3 import Connection
from typing import Iterable, NoReturn

from .cdna_config import CDNAConfig
from .codon_table import CodonTable
from .codon_table_builder import CodonTableBuilder
from .common_cli import finalise
from .db import get_db_conn, init_db
from .enums import SrcType
from .errors import SequenceNotFound
from .experiment_meta import ExperimentMeta
from .loaders.cdna_annot import load_annot
from .loaders.cdna_experiment_config import CDNAExperimentConfig
from .loaders.cdna_targeton_config import CDNATargetonConfig
from .loaders.fasta import load_from_multi_fasta
from .meta_table import MetaTable
from .oligo_generation_info import OligoGenerationInfo
from .oligo_seq import OligoSeq
from .queries import insert_annot_pattern_variants, insert_pattern_variants, clear_per_contig_tables, clear_per_targeton_tables, is_meta_table_empty
from .seq import Seq
from .strings.dna_str import DnaStr
from .strings.strand import Strand
from .transcript import Transcript
from .variant import Variant, VariantT


def exit_on_critical_exception(ex: Exception, msg: str) -> NoReturn:
    logging.critical(ex.args[0])
    logging.critical(msg)
    sys.exit(1)


def load_seqs(fp: str, ids: Iterable[str]) -> dict[str, Seq]:
    # TODO: use the actual targeton start position instead!
    return {
        seq_id: Seq(1, DnaStr(seq))
        for seq_id, seq in load_from_multi_fasta(fp, ids).items()
    }


def proc_targeton(
    conn: Connection,
    config: CDNAConfig,
    codon_table: CodonTable,
    seq: Seq,
    transcript: Transcript | None,
    targeton_cfg: CDNATargetonConfig
) -> OligoGenerationInfo:

    # Clear database tables
    clear_per_contig_tables.execute(conn)
    clear_per_targeton_tables.execute(conn)

    # Get targeton sequence
    targeton_seq = seq.subseq(targeton_cfg.ref, rel=False, prev_nt=True)

    # Get mutable region sequence
    r = targeton_cfg.region_2

    def get_vars_in_region(a: Iterable[VariantT]) -> list[VariantT]:
        return [x for x in a if x.pos in r and x.ref_end in r]

    r_seq = seq.subseq(r, rel=False)
    if transcript and not transcript.is_empty:
        # TODO: Load the single exon into the database for annotation purposes?
        cds = transcript.exons[0]
        if r.overlaps(cds):
            if r not in cds:
                raise ValueError(f"Partial CDS region {r}!")

            r_seq = transcript.get_cds_seq(seq, 0, r)

    # Generate pattern variants
    vars, annot_vars = targeton_cfg.mutator_collection.get_variants(
        codon_table, r_seq)

    vars = get_vars_in_region(vars)
    annot_vars = get_vars_in_region(annot_vars)

    def get_oligo(x: Variant) -> OligoSeq:
        # TODO: verify the need for a reference start!
        return OligoSeq.from_ref(targeton_seq, x, ref_start=x.pos)

    # Register the pattern variants into the database
    insert_pattern_variants(conn, list(map(get_oligo, vars)))
    insert_annot_pattern_variants(conn, list(map(get_oligo, annot_vars)))

    exp_meta = ExperimentMeta(
        config.species,
        config.assembly,
        config.min_length,
        config.max_length)

    targeton_name = f"{targeton_cfg.seq_id}_{targeton_cfg.get_hash()}"

    if not is_meta_table_empty(conn):
        MetaTable(
            SrcType.CDNA,
            targeton_cfg.seq_id,
            config,
            exp_meta,
            config.get_options(),
            None,
            '',
            Strand('+'),
            targeton_seq,
            targeton_seq,
            [],
            transcript.info if transcript else None).to_csv(conn, targeton_name)

    return OligoGenerationInfo()


def run_cdna(config: CDNAConfig) -> None:
    options = config.get_options()

    # Load targeton configurations
    try:
        targetons = CDNAExperimentConfig.load(config.oligo_info_fp)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA targeton file!")

    # Load codon table
    codon_table = CodonTableBuilder.load(
        config.get_codon_table_file_path()).build(Strand('+'))

    cds_info = {}
    try:
        # Load cDNA sequences
        seqs = load_seqs(config.ref_fasta_fp, targetons.sequence_ids)

        # Load transcript sequences
        if config.annot_fp:
            cds_info = load_annot(config.annot_fp, targetons.sequence_ids)

    except SequenceNotFound:
        sys.exit(1)
    except ValueError as ex:
        exit_on_critical_exception(ex, "Failed to load cDNA annotation file!")

    with get_db_conn() as conn:

        # Initialise in-memory database
        init_db(conn)

        proc_targeton_f = partial(
            proc_targeton,
            conn,
            config,
            codon_table)

        stats = OligoGenerationInfo()
        info: OligoGenerationInfo
        for targeton_cfg in targetons.targeton_configs:
            try:
                info = proc_targeton_f(
                    seqs[targeton_cfg.seq_id],
                    cds_info.get(targeton_cfg.seq_id),
                    targeton_cfg)

            except ValueError as ex:
                exit_on_critical_exception(ex, "Failed to generate oligonucleotides!")

            stats.update(info)

    finalise(config, stats)
