/*
================================================================================
    STAR - Spliced Transcripts Alignment to a Reference
================================================================================
*/

process STAR_ALIGN {

    tag "${meta.id}"
    label 'process_high'

    publishDir { "${params.outdir}/star/${meta.id}" },
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.bam')     ? "bam/$filename" :
            filename.endsWith('.bai')     ? "bam/$filename" :
            filename.endsWith('Log.final.out') ? "log/$filename" :
            filename.endsWith('.tab')     ? "spliceJunctions/$filename" : null
        }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*.sortedByCoord.out.bam"),       emit: bam
    tuple val(meta), path("*.sortedByCoord.out.bam.bai"),   emit: bai
    tuple val(meta), path("*Log.final.out"),                 emit: log_final
    tuple val(meta), path("*Log.out"),                       emit: log_out
    tuple val(meta), path("*SJ.out.tab"),                    emit: sj
    path  "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reads     = meta.single_end ? "${r1}" : "${r1} ${r2}"
    def two_pass  = params.star_two_pass ? "--twopassMode Basic" : ""
    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${params.star_index} \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFileNamePrefix ${prefix}. \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outBAMsortingThreadN ${task.cpus} \
        --limitBAMsortRAM 1500000000 \
        --runDirPerm All_RWX \
        ${two_pass} \
        ${args}


    # Index the BAM file
    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    # Rename for clarity
    mv ${prefix}.Aligned.sortedByCoord.out.bam \\
       ${prefix}.sortedByCoord.out.bam
    mv ${prefix}.Aligned.sortedByCoord.out.bam.bai \\
       ${prefix}.sortedByCoord.out.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
