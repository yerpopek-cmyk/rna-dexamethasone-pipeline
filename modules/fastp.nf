/*
================================================================================
    FASTP - Fast all-in-one FASTQ preprocessing
================================================================================
*/

process FASTP {

    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/fastp/${meta.id}" },
        mode: 'copy',
        pattern: "*.{html,json,log}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_R1.trimmed.fastq.gz"), path("*_R2.trimmed.fastq.gz"), emit: trimmed
    path  "*.html",                                                                  emit: html
    path  "*.json",                                                                  emit: json
    path  "*.log",                                                                   emit: log
    path  "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def paired  = meta.single_end ? "" : "-I ${r2} -O ${prefix}_R2.trimmed.fastq.gz"
    """
    fastp \\
        -i ${r1} \\
        -o ${prefix}_R1.trimmed.fastq.gz \\
        ${paired} \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --length_required ${params.fastp_min_length} \\
        --qualified_quality_phred ${params.fastp_quality} \\
        --html ${prefix}_fastp.html \\
        --json ${prefix}_fastp.json \\
        --correction \\
        --overrepresentation_analysis \\
        ${args} \\
        2> ${prefix}_fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | head -1 | sed 's/fastp //')
    END_VERSIONS
    """
}
