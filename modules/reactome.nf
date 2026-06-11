/*
================================================================================
    REACTOME_PA - Reactome pathway enrichment analysis (ORA)
================================================================================
*/

process REACTOME_ORA {

    label 'process_r'

    publishDir "${params.outdir}/pathway_analysis/reactome_ora",
        mode: 'copy'

    input:
    path deg_sig

    output:
    path "reactome_ORA_results.csv", emit: results
    path "ora_plots",                 emit: plots
    path "versions.yml",              emit: versions

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/bin/reactome.R \\
        --deg       ${deg_sig} \\
        --organism  ${params.organism} \\
        --pvalue    ${params.pvalue_cutoff} \\
        --lfc       ${params.lfc_threshold} \\
        --outdir    .

    mkdir -p ora_plots
    mv plots ora_plots/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //' | cut -d' ' -f1)
        bioconductor-reactomepa: \$(Rscript -e "cat(as.character(packageVersion('ReactomePA')))")
    END_VERSIONS
    """
}
