/*
================================================================================
    GSEA - Gene Set Enrichment Analysis via gsePathway (ReactomePA)
================================================================================
*/

process REACTOME_GSEA {

    label 'process_r'

    publishDir "${params.outdir}/pathway_analysis/reactome_gsea",
        mode: 'copy'

    input:
    path deg_all

    output:
    path "reactome_GSEA_results.csv", emit: results
    path "gsea_plots",                 emit: plots
    path "versions.yml",               emit: versions

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/bin/gsea.R \\
        --deg       ${deg_all} \\
        --organism  ${params.organism} \\
        --pvalue    ${params.pvalue_cutoff} \\
        --outdir    .

    mkdir -p gsea_plots
    mv plots gsea_plots/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //' | cut -d' ' -f1)
        bioconductor-reactomepa: \$(Rscript -e "cat(as.character(packageVersion('ReactomePA')))")
    END_VERSIONS
    """
}
