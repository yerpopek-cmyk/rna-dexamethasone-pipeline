/*
================================================================================
    REPORT - Render final RMarkdown HTML report
================================================================================
*/

process RENDER_REPORT {

    label 'process_r'

    publishDir "${params.outdir}/report",
        mode: 'copy'

    input:
    path deseq2_plots
    path ora_plots
    path gsea_plots

    output:
    path "*.html",        emit: report
    path "versions.yml",  emit: versions

    script:
    """
    Rscript ${projectDir}/bin/report.R \\
        --deseq2_dir   ${deseq2_plots} \\
        --ora_dir      ${ora_plots} \\
        --gsea_dir     ${gsea_plots} \\
        --rmd          ${projectDir}/report/RNAseq_report.Rmd \\
        --title        "${params.report_title}" \\
        --outdir       .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //' | cut -d' ' -f1)
        r-rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
    END_VERSIONS
    """
}
