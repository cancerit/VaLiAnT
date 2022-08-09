/* 
 * VaLiAnT module
 */

nextflow.enable.dsl = 2

// Defaults
params.valiant_config = ''

process valiant_sge {
    label "VaLiAnT"
    publishDir "${params.outdir}/valiant_sge", mode:'move'

    input:
    path(valiant_config)

    output:
    path("*_meta_excluded.csv"), emit: excluded_metadata_csv
    path("*_meta.csv"),          emit:metadata_csv
    path("*_unique.csv"),        emit:unique_csv
    path("*_pam.vcf"),           emit:pam_vcf
    path("*_ref.vcf"),           emit:ref_vcf
    path("config.json"),         emit:config_json
    path("ref_sequences.csv"),   emit:ref_sequences_csv

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