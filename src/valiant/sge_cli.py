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

from .common_cli import common_params, existing_file
from .sge_config import SGEConfig
from .sge_proc import run_sge


@click.command()
@common_params
@click.option('--gff', 'gff_fp', type=existing_file, help="Annotation GFF file path")
@click.option('--bg', 'bg_fp', type=existing_file, help="Background variant VCF file path")
@click.option('--bg-mask', 'mask_bg_fp', type=existing_file, help="Background variant excluded ranges BED file path")
@click.option('--pam', 'pam_fp', type=existing_file, help="PAM protection VCF file path")
@click.option('--vcf', 'vcf_fp', type=existing_file, help="Custom variant VCF manifest file path")
@click.option(
    '--sequences-only',
    is_flag=True,
    help="Fetch reference sequences and quit")
@click.option(
    '--revcomp-minus-strand',
    is_flag=True,
    help="Include reverse complement in oligonucleotide if reference is on minus strand")
@click.option(
    '--force-bg-ns',
    'force_bg_ns',
    is_flag=True,
    help="Allow non-synonymous background variants")
@click.option(
    '--force-bg-indels',
    'force_bg_fs',
    is_flag=True,
    help="Allow frame-shifting background variants")
def sge(

    # Input files
    oligo_info_fp: str,
    ref_fasta_fp: str,
    codon_table_fp: str | None,
    gff_fp: str | None,
    bg_fp: str | None,
    pam_fp: str | None,
    vcf_fp: str | None,

    # Output directory
    output_dir: str,

    # Metadata
    species: str,
    assembly: str,

    # Adaptor sequences
    adaptor_5: str | None,
    adaptor_3: str | None,

    # Background variants
    force_bg_ns: bool,
    force_bg_fs: bool,
    mask_bg_fp: str | None,

    # Actions
    sequences_only: bool,
    revcomp_minus_strand: bool,
    max_length: int,
    min_length: int

) -> None:
    """
    SGE oligonucleotide generation tool

    \b
    OLIGO_INFO is the BED-like input file path
    REF_FASTA is the reference genome FASTA file path
    OUTPUT is the output directory path
    SPECIES will be included in the metadata
    ASSEMBLY will be included in the metadata
    """

    run_sge(
        SGEConfig(
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
            revcomp_minus_strand=revcomp_minus_strand,
            force_bg_ns=force_bg_ns,
            force_bg_fs=force_bg_fs,
            gff_fp=gff_fp,
            bg_fp=bg_fp,
            pam_fp=pam_fp,
            vcf_fp=vcf_fp,
            mask_bg_fp=mask_bg_fp),
        sequences_only)
