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
include { FASTQC            } from '../../modules/nf-core/fastqc/main'
include { MULTIQC           } from '../../modules/nf-core/multiqc/main'

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
include { FASTA_FASTQ_PROTEIN_ALIGNMENT_DIAMOND } from '../../subworkflows/local/fasta_fastq_protein_alignment_diamond/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAGANALYZER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    //
    // extract relevant files from samplesheet
    //
    ch_fastq    = ch_samplesheet.map { it[0], it[1] }
    ch_fna      = ch_samplesheet.map { it[0], it[2] }
    ch_faa      = ch_samplesheet.map { it[0], it[3] }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


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
    //
    // SUBWORKFLOW: Predict genes for input nucleotide sequences
    //



    /*
    ------------------------------------------------------------------------------
        Cluster and dereplicate gene sequences
    ------------------------------------------------------------------------------
    */
    //
    // SUBWORKFLOW: Predict genes for input nucleotide sequences
    //


    /*
    ------------------------------------------------------------------------------
        Read alignment
    ------------------------------------------------------------------------------
    */
    //
    // SUBWORKFLOW: Align reads to protein catalog
    //
    ch_fastq_diamond_results = FASTA_FASTQ_PROTEIN_ALIGNMENT_DIAMOND (
        ch_fastq,
        ch_fna,
        "txt",
        "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    ).diamond_results_txt
    ch_versions = ch_versions.mix(FASTA_FASTQ_PROTEIN_ALIGNMENT_DIAMOND.out.versions())

    //
    // MODULE: Assign multi-mapping reads to by likelihood inference
    //

    //
    // MODULE: Create a file containing the abundance of each gene in each sample
    //


    /*
    ------------------------------------------------------------------------------
        Create co-abundant gene groups (CAGs) using alignment information
    ------------------------------------------------------------------------------
    */
    //
    // SUBWORKFLOW: Create and refine CAGs
    //

    //
    // MODULE: Create a file containing the abundance of each CAG in each sample
    //


    /*
    ------------------------------------------------------------------------------
        Annotate CAGs using eggnog mapper
    ------------------------------------------------------------------------------
    */






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
