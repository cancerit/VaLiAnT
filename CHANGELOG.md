# Changelog

## [UNRELEASED]

### Added

- Support for background variant
- Support for out-of-frame CDS targeton regions whose 5' and 3' extensions span both adjacent and distal bases

### Fixed

- Fix codon tables loaded from files where the codons are not grouped by amino acid

## [3.0.1] - 2023-07-18

### Fixed

- Fix custom variants being ignored when no targeton region is set as mutable
- Fix error on mutations that end but do not start in codons affected by PAM protection edits
- Metadata table: fix MAVE-HGVS strings (`mave_nt` field) to correctly report as deletion-insertions those deletions that occur in codons affected by PAM protection edits where the edit itself is preserved

## [3.0.0] - 2022-10-11

### Added

- CLI: add `--config` option to main command
- CLI: add `--min-length` option to both subcommands
- Metadata table: add `mseq_no_adapt` field
- Metadata table: add `pam_mut_annot` field
- Metadata table: add `pam_mut_sgrna_id` field
- Metadata table: add `mave_nt` field
- Metadata table: add `mave_nt_ref` field
- Metadata table: add `vcf_var_in_const` field
- Output files: generate additional VCF file where `REF` does not include PAM variants
- Output files: generate JSON file reporting the execution parameters

### Changed

- CLI: enable main command to run from configuration
- Output files: extend the `REF` of the entries in the output VCF to include the PAM variants occurring at the liminal codons, if any

### Fixed

- Fix the `REF` of custom substitutions and deletions to include PAM variants, if any (this may affect oligonucleotide names)

## [2.0.0] - 2021-07-12

### Added

- Add cDNA DMS library design as `cdna` subcommand

### Changed

- Move SGE library design to `sge` subcommand
- Lock dependency versions
- Make unique oligonucleotide name selection order-invariant

### Fixed

- Fix discarded metadata column order and filtering

## [1.0.0] - 2020-12-21

Initial release.
