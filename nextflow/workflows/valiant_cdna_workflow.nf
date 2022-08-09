/* 
 * Workflow to run VaLiAnT
 */

nextflow.enable.dsl = 2

// Default parameters
params.help = false
params.outdir = "."
params.valiant_config=""

// Modules
include { valiant_cdna } from '../modules/run_valiant_cdna.nf'

// Define usage
def helpMessage() {
  log.info """
        VaLiAnT workflow
        -------------------------------------------------------
        Usage:
        The typical command for running the workflow is as follows:
        nextflow -C ${PWD}/config/nextflow.config run ${PWD}/workflows/valiant_workflow.nf --valiant_config <path>

        Mandatory arguments:
          --valiant_config PATH  Path for VaLiAnT JSON-formatted config file.
        Optional arguments:
          --help                         This usage statement.
        """
}

// Print usage 
if (params.help) {
  helpMessage()
  exit 0
}

// Output directory validation
// If none is set, will default to current working directory
if (params.outdir) {
  if (!file(params.outdir, type: 'dir', checkIfExists: true)) {
  exit 1, "Output directory does not exist: ${params.outdir}."
  }
}

// VaLiAnT config validation
if (!params.valiant_config) {
  exit 1, "Undefined --valiant_config parameter. Please provide path to VaLiAnT JSON-formatted config file."
} 
if (!file(params.valiant_config, checkIfExists: true)) {
  exit 1, "VaLiAnT JSON-formatted config file does not exist: ${params.valiant_config}."
}

// Logging
log.info 'Starting workflow.....'
log.info "VaLiAnT config: ${params.valiant_config}"

workflow {
  valiant_cdna(params.valiant_config)
}