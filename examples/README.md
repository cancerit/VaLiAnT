# Example libraries

For more information on these examples, please see the [VaLiAnT Wiki](https://github.com/cancerit/VaLiAnT/wiki).

These instructions assume you have installed the `valiant` package in your local environment. The scripts can be easily adapted to use the Docker or Singularity image by mounting the necessary paths as described in the [main documentation](../README.md).

Once the examples are run, their outputs can be compared with the expected results via the corresponding `check_*.sh` scripts. The presence or absence of certain types of output file will depend on the execution parameters. The validation scripts depend on the `md5sum` utility (with `md5` as a fallback).

On macOS, the `md5sum` tool can be installed, *e.g.*, via [Homebrew](https://brew.sh/):

```sh
brew install coreutils
```

## SGE

To run the examples, first unpack the reference sequence:

```sh
cd sge
./unpack_reference.sh
```

### BRCA1

To generate the example BRCA1 libraries (the corresponding output directories will be created if they don't exist):

```sh
cd brca1
./run_brca1_nuc.sh  # -> brca1_nuc_output/
./run_brca1_pep.sh  # -> brca1_pep_output/
./check_brca1_nuc.sh
./check_brca1_pep.sh
```

The expected outputs (the order of the oligonucleotides may differ) are in the `brca1_nuc_output_exp` and `brca1_pep_output_exp` directories.

### BRCA1 prime editing

To run the examples, first unpack the input files:

```sh
cd brca1_prime_editing
./unpack_inputs.sh
```

To generate the example BRCA1 libraries (the corresponding output directories will be created if they don't exist):

```sh
./run_prime_a.sh  # -> output_a/
./run_prime_b.sh  # -> output_b/
./check_prime_a.sh
./check_prime_b.sh
```

The expected outputs (the order of the oligonucleotides may differ) are in the `output_a_exp` and `output_b_exp` directories.

### NXRL

This example applies background variants.

```sh
cd nxrl
./run.sh  # -> output/
./check.sh
```

The expected outputs (the order of the oligonucleotides may differ) are in the `output_exp` directory.

### cDNA DMS

```sh
cd cdna
./run.sh  # -> output/
./check.sh
```

The expected outputs (the order of the oligonucleotides may differ) are in the `output_exp` directory.
