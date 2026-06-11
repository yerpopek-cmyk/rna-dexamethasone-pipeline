/*
================================================================================
    DIFFERENTIAL_EXPRESSION - Subworkflow: DESeq2 + pathway analysis + report
================================================================================
*/

include { DESEQ2         } from '../modules/deseq2'
include { REACTOME_ORA   } from '../modules/reactome'
include { REACTOME_GSEA  } from '../modules/gsea'
include { RENDER_REPORT  } from '../modules/report'

workflow DIFF_EXPRESSION {

    take:
    counts    // path: gene_counts.txt
    metadata  // path: metadata.csv

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: DESeq2 differential expression
    //
    DESEQ2(counts, metadata)
    ch_versions = ch_versions.mix(DESEQ2.out.versions)

    //
    // MODULE: Reactome Over-Representation Analysis (ORA)
    //
    REACTOME_ORA(DESEQ2.out.deg_sig)
    ch_versions = ch_versions.mix(REACTOME_ORA.out.versions)

    //
    // MODULE: Reactome Gene Set Enrichment Analysis (GSEA)
    //
    REACTOME_GSEA(DESEQ2.out.deg_table)
    ch_versions = ch_versions.mix(REACTOME_GSEA.out.versions)

    //
    // MODULE: Render final HTML report
    //
    RENDER_REPORT(
        DESEQ2.out.plots,
        REACTOME_ORA.out.plots,
        REACTOME_GSEA.out.plots
    )
    ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)

    emit:
    deg_table     = DESEQ2.out.deg_table          // path: DEG_results.csv
    deg_sig       = DESEQ2.out.deg_sig            // path: DEG_significant.csv
    norm_counts   = DESEQ2.out.norm_counts        // path: normalized_counts.csv
    pathway_table = REACTOME_ORA.out.results      // path: reactome_ORA_results.csv
    gsea_table    = REACTOME_GSEA.out.results     // path: reactome_GSEA_results.csv
    html_report   = RENDER_REPORT.out.report      // path: *.html
    versions      = ch_versions
}
