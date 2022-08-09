# Nextflow module and workflow for VaLiAnT

## Installing Nextflow

Set the version of Nextflow to install:

```
nextflow_version='22.04.3'
```

Download Nextflow with all dependencies:

```
wget https://github.com/nextflow-io/nextflow/releases/download/v${nextflow_version}/nextflow-${nextflow_version}-all
mv nextflow-${nextflow_version}-all nextflow
chmod +x nextflow
```

### Binding local paths

It may be necessary to bind or mount your local directories so that Docker or Singularity can see the input data and output directory. Below are some examples which can be modified and added to the Nextflow config.

For Docker:
```
docker {
    enabled = true
    runOptions = '-v /your/local/path:/container/path'
}
```

For Singularity:

```
singularity {
    enabled = true
    runOptions = '--bind /your/local/path'
}
```

Alternatively, the paths can be set using an environment variable:

```
env {
  VALIANT_REPO = "/your/local/path:/container/path"
}

docker {
  enabled = true
  runOptions = "-v $env.VALIANT_REPO"
}
```

## Using Nextflow to run example datasets with VaLiAnT

Before running any commands, you will need to make sure:

* **[Docker](https://www.docker.com/) is installed**
* You have `nextflow` available in your PATH (e.g. `export PATH=$PATH:/path/to/directory/with/nextflow`)
* You have updated the `VALIANT_REPO` variable in `nextflow/examples/example.config` to reflect the full path to the VaLiAnT repository on your machine.

All outputs will be written the directory from which the pipeline is run unless `--outdir` is specified.

For the following examples, we'll set an environment variable to make the commands agnostic.

```
LOCAL_VALIANT_REPO=/path/to/valiant/repo
```

### SGE

Command to run:

```
nextflow -C ${LOCAL_VALIANT_REPO}/nextflow/examples/example.config run ${LOCAL_VALIANT_REPO}/nextflow/workflows/valiant_sge_workflow.nf --valiant_config ${LOCAL_VALIANT_REPO}/nextflow/examples/sge/valiant.json 
```

Results will be written to a subdirectory called `valiant_sge` in the current directory.

### cDNA

Command to run:

```
nextflow -C ${LOCAL_VALIANT_REPO}/nextflow/examples/example.config run ${LOCAL_VALIANT_REPO}/nextflow/workflows/valiant_cdna_workflow.nf --valiant_config ${LOCAL_VALIANT_REPO}/nextflow/examples/cdna/valiant.json 
```

Results will be written to a subdirectory called `valiant_cdna` in the current directory.