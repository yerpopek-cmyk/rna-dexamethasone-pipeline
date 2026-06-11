#!/usr/bin/env Rscript
# =============================================================================
# Render Final RMarkdown HTML Report
# Project: RNA-seq Dexamethasone Effects (GSE52778)
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(rmarkdown)
})

option_list <- list(
    make_option("--deseq2_dir", type="character", help="DESeq2 results directory"),
    make_option("--ora_dir",    type="character", help="Reactome ORA results directory"),
    make_option("--gsea_dir",   type="character", help="Reactome GSEA results directory"),
    make_option("--rmd",        type="character", help="Path to .Rmd template"),
    make_option("--title",      type="character", default="RNA-seq Analysis Report",
                                help="Report title"),
    make_option("--outdir",     type="character", default=".", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

message("=== Rendering HTML Report ===")
message("  Rmd template : ", opt$rmd)
message("  DESeq2 dir   : ", opt$deseq2_dir)
message("  ORA dir      : ", opt$ora_dir)
message("  GSEA dir     : ", opt$gsea_dir)
message("  Output dir   : ", opt$outdir)

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

output_file <- file.path(normalizePath(opt$outdir), "RNAseq_analysis_report.html")

rmarkdown::render(
    input          = opt$rmd,
    output_file    = output_file,
    output_format  = "html_document",
    params         = list(
        deseq2_dir = normalizePath(opt$deseq2_dir),
        ora_dir    = normalizePath(opt$ora_dir),
        gsea_dir   = normalizePath(opt$gsea_dir),
        title      = opt$title
    ),
    envir   = new.env(parent=globalenv()),
    quiet   = FALSE
)

message("\nReport rendered: ", output_file)
