/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/FUNCTIONS/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// FUNCTIONS: Imported from nf-core plugins
//
include { paramsSummaryMap       } from 'plugin/nf-validation'

//
// MODULES: Imported from nf-core/modules
//
include { FASTQC                } from '../../modules/nf-core/fastqc/main'
include { PRODIGAL              } from '../../modules/nf-core/prodigal/main'
include { METAEUK_EASYPREDICT   } from '../../modules/nf-core/metaeuk/easypredict/main'
include { MULTIQC               } from '../../modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS: Imported from nf-core/modules
//
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_caganalyzer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULES: Local modules
//
include { PRODIGALGV    } from '../../modules/nf-core/prodigalgv/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAGANALYZER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // extract relevant files from samplesheet
    //
    ch_fastq    = ch_samplesheet.map { it[0], it[1] }
    ch_fna      = ch_samplesheet.map { it[0], it[2] }
    ch_faa      = ch_samplesheet.map { it[0], it[3] }


    /*
    ------------------------------------------------------------------------------
        Read quality assessment
    ------------------------------------------------------------------------------
    */
    //
    // MODULE: Quality assess input reads
    //
    FASTQC (
        ch_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    /*
    ------------------------------------------------------------------------------
        OPTIONAL: Gene prediction
    ------------------------------------------------------------------------------
    */
    // if params.gene_prediction == "prodigal", use prodigal
    if ( params.gene_prediction_tool == "prodigal" ) {
        //
        // MODULE: Predict genes for input nucleotide sequences using Prodigal (default)
        //
        ch_gene_predictions_gff = PRODIGAL ( ch_fna, "gff" ).out.gene_annotations
        ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())
    } else if ( params.gene_prediction_tool == "prodigal-gv" ) {
        //
        // MODULE: Predict genes for input nucleotide sequences using Prodigal-gv (for viruses/phages)
        //
        ch_gene_predictions_gff = PRODIGALGV ( ch_fna, "gff" ).out.gene_annotations
        ch_versions = ch_versions.mix(PRODIGALGV.out.versions.first())
    } else if ( params.gene_prediction_tool == "metaeuk") {
        //
        // MODULE: Predict genes for input nucleotide sequences using metaEuk (for eukaryotes)
        //
        ch_gene_predictions_gff = METAEUK_EASYPREDICT ( ch_fna, "gff" ).out.gff
        ch_versions = ch_versions.mix(METAEUK_EASYPREDICT.out.versions.first())
    }






    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
