/*
================================================================================
    DESEQ2 - Differential expression analysis
================================================================================
*/

process DESEQ2 {

    label 'process_r'

    publishDir "${params.outdir}/deseq2",
        mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "DEG_results.csv",         emit: deg_table
    path "DEG_significant.csv",     emit: deg_sig
    path "normalized_counts.csv",   emit: norm_counts
    path "vst_counts.csv",          emit: vst_counts
    path "deseq2_plots",           emit: plots
    path "rds",                    emit: rds_objects
    path "versions.yml",           emit: versions

    script:
    """
    mkdir -p plots rds

    Rscript ${projectDir}/bin/deseq2.R \\
        --counts      ${counts} \\
        --metadata    ${metadata} \\
        --pvalue      ${params.pvalue_cutoff} \\
        --lfc         ${params.lfc_threshold} \\
        --min_count   ${params.deseq2_min_count} \\
        --lfc_shrink  ${params.deseq2_lfc_shrink} \\
        --outdir      .

    mkdir -p deseq2_plots
    mv plots deseq2_plots/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //' | cut -d' ' -f1)
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
