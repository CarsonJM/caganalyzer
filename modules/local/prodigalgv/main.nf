process PRODIGALGV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prodigal-gv:2.11.0--he4a0461_2':
        'biocontainers/prodigal-gv:2.11.0--he4a0461_2' }"

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}.gz"),    emit: gene_annotations
    tuple val(meta), path("${prefix}.fna.gz"),                 emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa.gz"),                 emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt.gz"),             emit: all_gene_annotations
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${genome} | prodigal-gv \\
        $args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    pigz -nm ${prefix}.fna
    pigz -nm ${prefix}.${output_format}
    pigz -nm ${prefix}.faa
    pigz -nm ${prefix}_all.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal-gv: \$(prodigal-gv -v 2>&1 | sed -n 's/Prodigal-gv V\\(.*\\):.*/\\1/p')
        pigz: \$(pigz -V 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fna.gz
    touch ${prefix}.${output_format}.gz
    touch ${prefix}.faa.gz
    touch ${prefix}_all.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal-gv: \$(prodigal-gv -v 2>&1 | sed -n 's/Prodigal-gv V\\(.*\\):.*/\\1/p')
        pigz: \$(pigz -V 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

}
