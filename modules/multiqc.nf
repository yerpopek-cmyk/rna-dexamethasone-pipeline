/*
================================================================================
    MULTIQC - Aggregate bioinformatics results across many samples
================================================================================
*/

process MULTIQC {

    label 'process_low'

    publishDir "${params.outdir}/multiqc",
        mode: 'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "*_plots",              emit: plots,   optional: true
    path "versions.yml",         emit: versions

    script:
    def args     = task.ext.args   ?: ''
    def config   = file("${projectDir}/assets/multiqc_config.yml").exists()
                   ? "--config ${projectDir}/assets/multiqc_config.yml" : ''
    def logo     = file("${projectDir}/assets/multiqc_logo.png").exists()
                   ? "--cl-config 'custom_logo: ${projectDir}/assets/multiqc_logo.png'" : ''
    """
    mkdir -p /tmp/multiqc_out
    multiqc \\
        --force \\
        --title "RNA-seq Dexamethasone Analysis" \\
        --comment "Quality control report for GSE52778 dataset" \\
        ${config} \\
        ${logo} \\
        ${args} \\
        -o /tmp/multiqc_out \\
        .

    cp -r /tmp/multiqc_out/* .
    rm -rf /tmp/multiqc_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
