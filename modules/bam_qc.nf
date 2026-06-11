/*
================================================================================
    BAM QC - SAMtools flagstat + stats + idxstats in one process to minimize file copies
================================================================================
*/

process BAM_QC {

    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/bam_qc/${meta.id}" },
        mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path  "versions.yml",                 emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools flagstat \\
        --threads ${task.cpus} \\
        ${bam} > ${prefix}.flagstat

    samtools stats \\
        --threads ${task.cpus} \\
        ${bam} > ${prefix}.stats

    samtools idxstats ${bam} > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
