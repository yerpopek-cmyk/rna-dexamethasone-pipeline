/*
================================================================================
    FASTQC - Sequence quality assessment
================================================================================
*/

process FASTQC {

    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/fastqc/${meta.id}" },
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.zip')  ? "zip/$filename" :
            filename.endsWith('.html') ? "html/$filename" : null
        }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def reads    = meta.single_end ? "${r1}" : "${r1} ${r2}"
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC //')
    END_VERSIONS
    """
}
