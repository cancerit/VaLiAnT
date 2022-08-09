/* 
 * VaLiAnT module
 */

nextflow.enable.dsl = 2

// Defaults
params.valiant_config = ''

process valiant_cdna {
    label "VaLiAnT"
    publishDir "${params.outdir}/valiant_cdna", mode:'move'

    input:
    path(valiant_config)

    output:
    path("*_meta.csv"),          emit:metadata_csv
    path("*_unique.csv"),        emit:unique_csv
    path("config.json"),         emit:config_json

    script:
    if( !valiant_config.exists() ) {
        exit 1, "VaLiAnT config does not exist: ${valiant_config}"
    }
    else {
        """
        valiant -c ${valiant_config}
        """	
    }
}