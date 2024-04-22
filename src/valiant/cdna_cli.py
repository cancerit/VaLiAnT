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

import click

from .cdna_config import CDNAConfig
from .cdna_proc import run_cdna
from .common_cli import common_params, existing_file


@click.command()
@common_params
@click.option('--annot', 'annot_fp', type=existing_file, help="cDNA annotation file path")
def cdna(

    # Input files
    oligo_info_fp: str,
    ref_fasta_fp: str,
    codon_table_fp: str | None,
    annot_fp: str | None,

    # Output directory
    output_dir: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: str | None,
    adaptor_3: str | None,

    # Actions
    max_length: int,
    min_length: int

):
    """
    cDNA DMS oligonucleotide generation tool

    \b
    OLIGO_INFO is the BED-like input file path
    REF_FASTA is the cDNA multi-FASTA file path
    OUTPUT is the output directory path
    SPECIES will be included in the metadata
    ASSEMBLY will be included in the metadata
    """

    run_cdna(CDNAConfig(
        species=species,
        assembly=assembly,
        adaptor_5=adaptor_5,
        adaptor_3=adaptor_3,
        min_length=min_length,
        max_length=max_length,
        codon_table_fp=codon_table_fp,
        oligo_info_fp=oligo_info_fp,
        ref_fasta_fp=ref_fasta_fp,
        output_dir=output_dir,
        annot_fp=annot_fp))
