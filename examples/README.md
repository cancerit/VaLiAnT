# Example libraries

These instructions assume you have installed the Valiant package in your local environment. The scripts can be easily adapted to use the Docker or Singularity image by mounting the necessary paths as described in the [main documentation](../README.md).

### SGE

To run the examples, first unpack the reference sequence:

```sh
cd sge
./unpack_reference.sh
```

To generate the example BRCA1 libraries (the corresponding output directories will be created if they don't exist):

```sh
./run_brca1_nuc.sh  # -> brca1_nuc_output/
./run_brca1_pep.sh  # -> brca1_pep_output/
```

The expected outputs (the order of the oligonucleotides may differ) are in the `brca1_nuc_output_exp` and `brca1_pep_output_exp` directories.

### cDNA DMS

```sh
cd cdna
./run_brca1_cdna.sh  # -> brca1_cdna_output/
```

The expected outputs (the order of the oligonucleotides may differ) are in the `brca1_cdna_output_exp` directory.
