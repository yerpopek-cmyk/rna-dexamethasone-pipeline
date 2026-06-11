/*
================================================================================
    FEATURECOUNTS - Counting reads over genomic features
================================================================================
*/

process FEATURECOUNTS {

    label 'process_medium'

    publishDir "${params.outdir}/featurecounts",
        mode: 'copy'

    input:
    val bam_files
    val gtf

    output:
    path "gene_counts.txt",         emit: counts
    path "gene_counts.txt.summary", emit: summary
    path "versions.yml",             emit: versions

    script:
    def paired_flag = "-p --countReadPairs"
    def bam_list = bam_files.join(' ')
    """
    # NOTE: This process runs outside Docker (container = null) to use a
    # host-compiled featureCounts binary. On WSL2, Docker-compiled binaries
    # can segfault due to kernel version mismatches. Compile from source:
    #   wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz
    #   tar xzf subread-2.0.6-source.tar.gz && cd subread-2.0.6-source/src
    #   make -f Makefile.Linux
    # Then add the bin/ directory to your PATH, or update the path below.
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o gene_counts.txt \\
        -t ${params.fc_feature_type} \\
        -g ${params.fc_group_attribute} \\
        -B \\
        -C \\
        --minOverlap 10 \\
        --fracOverlap 0.0 \\
        ${paired_flag} \\
        ${bam_list} \\
        2>&1 | tee featurecounts.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oP 'v\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}
