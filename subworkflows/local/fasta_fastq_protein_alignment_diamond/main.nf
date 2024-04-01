//
// Create and align to a DIAMOND database
//
include { DIAMOND_MAKEDB    } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX    } from '../../modules/nf-core/diamond/blastx/main'

workflow FASTA_FASTQ_PROTEIN_ALIGNMENT_DIAMOND {

    take:
    ch_fastq                    // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    ch_protein_fasta            // channel: [ val(meta), [ protein_fasta ] ]
    ch_diamond_output_format    // channel: val (string)
    ch_diamond_output_columns   // channel: val [qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Create DIAMOND database from protein fasta
    //
    ch_diamond_db   = DIAMOND_MAKEDB ( ch_protein_fasta, [], [], [] ).db
    ch_versions     = ch_versions.mix(DIAMOND_MAKEDB.out.versions.first())

    //
    // MODULE: Create DIAMOND database from protein fasta
    //
    ch_diamond_results  = DIAMOND_BLASTX ( ch_fastq, ch_diamond_db, ch_diamond_output_format, ch_diamond_output_columns ).out
    ch_versions         = ch_versions.mix(DIAMOND_MAKEDB.out.versions.first())

    emit:
    diamond_output  = ch_diamond_results    // channel: [ all_diamond_outputs ]
    versions        = ch_versions           // channel: [ versions.yml ]
}

