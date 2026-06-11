#!/usr/bin/env Rscript
# =============================================================================
# Gene Set Enrichment Analysis (GSEA) via ReactomePA::gsePathway
# Project: RNA-seq Dexamethasone Effects (GSE52778)
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ReactomePA)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(RColorBrewer)
})

# ── Argument parsing ────────────────────────────────────────────────────────

option_list <- list(
    make_option("--deg",      type="character", help="Full DEG results CSV (all genes, from deseq2.R)"),
    make_option("--organism", type="character", default="human", help="Organism for ReactomePA"),
    make_option("--pvalue",   type="double",    default=0.05,    help="p-value cutoff for GSEA"),
    make_option("--outdir",   type="character", default=".",     help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

set.seed(42)

message("=== Reactome GSEA Analysis ===")
message("Parameters:")
message("  DEG input : ", opt$deg)
message("  Organism  : ", opt$organism)
message("  p-value   : ", opt$pvalue)

plot_dir <- file.path(opt$outdir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

save_plot <- function(plot_obj, filename, width=12, height=10) {
    base <- file.path(plot_dir, filename)
    tryCatch({
        ggsave(paste0(base, ".pdf"), plot_obj, width=width, height=height)
        ggsave(paste0(base, ".png"), plot_obj, width=width, height=height, dpi=300)
    }, error=function(e) message("  Warning: could not save ", filename, ": ", e$message))
}

# ── Load all DEGs ────────────────────────────────────────────────────────────

message("\n[1/4] Loading all DEG results for ranked list...")
deg_all <- read.csv(opt$deg)
if (nrow(deg_all) > 0 && "gene_id" %in% colnames(deg_all)) {
    deg_all$gene_id <- sub("\\..*$", "", deg_all$gene_id)
}

# Remove NA log2FoldChange
deg_all <- deg_all %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    arrange(desc(log2FoldChange))

message("  Genes loaded: ", nrow(deg_all))

# ── ID Conversion ─────────────────────────────────────────────────────────────

message("\n[2/4] Converting gene IDs to Entrez for GSEA...")

gene_ids  <- deg_all$gene_id
is_ensembl <- grepl("^ENSG", gene_ids[1])
from_type  <- if(is_ensembl) "ENSEMBL" else "SYMBOL"

id_map <- bitr(
    gene_ids,
    fromType = from_type,
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db,
    drop     = TRUE
)

# Build ranked gene list: log2FC ranked descending
# Merge DEG results with Entrez IDs
deg_mapped <- deg_all %>%
    inner_join(id_map, by=setNames(from_type, "gene_id"))

# Handle duplicate Entrez IDs: keep entry with highest |log2FC|
deg_mapped <- deg_mapped %>%
    group_by(ENTREZID) %>%
    slice_max(order_by=abs(log2FoldChange), n=1, with_ties=FALSE) %>%
    ungroup()

# Create named numeric vector: Entrez ID -> log2FC
gene_list <- setNames(deg_mapped$log2FoldChange, deg_mapped$ENTREZID)
gene_list <- sort(gene_list, decreasing=TRUE)

message("  Mapped genes for GSEA: ", length(gene_list))
message("  Score range          : [", round(min(gene_list),2), ", ", round(max(gene_list),2), "]")

# ── Run GSEA ──────────────────────────────────────────────────────────────────

message("\n[3/4] Running Reactome GSEA (gsePathway)...")

gsea_result <- tryCatch(
    gsePathway(
        geneList      = gene_list,
        organism      = opt$organism,
        exponent      = 1,
        minGSSize     = 10,
        maxGSSize     = 500,
        eps           = 1e-10,
        pvalueCutoff  = opt$pvalue,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        seed          = 42,
        by            = "fgsea"
    ),
    error = function(e) {
        message("ERROR in gsePathway: ", e$message)
        return(NULL)
    }
)

if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    message("WARNING: No significant GSEA pathways found at p.adj < ", opt$pvalue)
    message("  Trying with relaxed cutoff (p.adj < 0.1)...")

    gsea_result <- tryCatch(
        gsePathway(
            geneList      = gene_list,
            organism      = opt$organism,
            pvalueCutoff  = 0.1,
            pAdjustMethod = "BH",
            verbose       = FALSE,
            seed          = 42
        ),
        error = function(e) { message("  Still failed: ", e$message); NULL }
    )
}

if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
    n_activated  <- sum(gsea_result@result$NES > 0)
    n_suppressed <- sum(gsea_result@result$NES < 0)
    message("  Significant pathways: ", nrow(gsea_result))
    message("    Activated  (NES > 0): ", n_activated)
    message("    Suppressed (NES < 0): ", n_suppressed)
} else {
    message("  No significant GSEA pathways found.")
}

# ── Export results ────────────────────────────────────────────────────────────

if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
    gsea_df <- as.data.frame(gsea_result)
    write.csv(gsea_df,
              file.path(opt$outdir, "reactome_GSEA_results.csv"),
              row.names=FALSE, quote=FALSE)
    message("  Saved: reactome_GSEA_results.csv (", nrow(gsea_df), " pathways)")
} else {
    write.csv(data.frame(),
              file.path(opt$outdir, "reactome_GSEA_results.csv"),
              row.names=FALSE)
}

# ── Plots ─────────────────────────────────────────────────────────────────────

message("\n[4/4] Generating GSEA plots...")

if (!is.null(gsea_result) && nrow(gsea_result) > 0) {

    n_show <- min(20, nrow(gsea_result))

    # 1. Ridgeplot (distribution of enrichment scores)
    message("  - Ridgeplot")
    tryCatch({
        p_ridge <- ridgeplot(gsea_result, showCategory=n_show, fill="p.adjust") +
            scale_fill_viridis_c(direction=-1, name="p.adjust") +
            labs(
                title    = "GSEA: Enrichment Score Distributions",
                subtitle = "Reactome pathway gene sets — ranked by NES",
                x        = "Enrichment score (log2FC distribution)",
                y        = "Pathway"
            ) +
            theme_bw(base_size=11) +
            theme(
                plot.title  = element_text(face="bold", size=14),
                axis.text.y = element_text(size=8)
            )
        save_plot(p_ridge, "GSEA_ridgeplot", width=14, height=max(8, n_show*0.5))
    }, error=function(e) message("  Warning: ridgeplot failed: ", e$message))

    # 2. Dotplot GSEA
    message("  - Dotplot")
    tryCatch({
        p_dot <- dotplot(gsea_result, showCategory=n_show, split=".sign") +
            facet_grid(. ~ .sign) +
            labs(title="GSEA: Enriched Pathways",
                 subtitle="Separated by activation/suppression") +
            theme_bw(base_size=11) +
            theme(
                plot.title  = element_text(face="bold", size=14),
                axis.text.y = element_text(size=8)
            )
        save_plot(p_dot, "GSEA_dotplot", width=16, height=10)
    }, error=function(e) message("  Warning: GSEA dotplot failed: ", e$message))

    # 3. NES barplot (top activated and suppressed)
    message("  - NES barplot")
    tryCatch({
        gsea_df2 <- as.data.frame(gsea_result) %>%
            arrange(desc(NES)) %>%
            mutate(
                Direction = ifelse(NES > 0, "Activated", "Suppressed"),
                Description = stringr::str_wrap(Description, 40)
            )

        top_n_each <- 10
        plot_df <- bind_rows(
            head(gsea_df2 %>% filter(Direction=="Activated"),  top_n_each),
            tail(gsea_df2 %>% filter(Direction=="Suppressed"), top_n_each)
        ) %>%
            mutate(Description = factor(Description, levels=rev(Description)))

        p_nes <- ggplot(plot_df,
                        aes(x=NES, y=Description, fill=Direction)) +
            geom_col(width=0.7, alpha=0.85) +
            geom_vline(xintercept=0, linewidth=0.5) +
            scale_fill_manual(values=c("Activated"="#D73027", "Suppressed"="#4575B4")) +
            labs(
                title    = "GSEA: Normalized Enrichment Scores",
                subtitle = paste0("Top ", top_n_each, " activated and suppressed Reactome pathways"),
                x        = "Normalized Enrichment Score (NES)",
                y        = NULL,
                fill     = NULL
            ) +
            theme_bw(base_size=11) +
            theme(
                plot.title   = element_text(face="bold", size=14),
                axis.text.y  = element_text(size=8),
                legend.position = "top"
            )
        save_plot(p_nes, "GSEA_NES_barplot", width=14, height=10)
    }, error=function(e) message("  Warning: NES barplot failed: ", e$message))

    # 4. Enrichment plots for top 4 pathways
    message("  - Enrichment running score plots (top 4)")
    tryCatch({
        top_ids <- head(gsea_result@result %>%
                            arrange(p.adjust) %>%
                            pull(ID), 4)

        for (i in seq_along(top_ids)) {
            path_id   <- top_ids[i]
            path_desc <- gsea_result@result$Description[
                gsea_result@result$ID == path_id][1]
            path_desc <- stringr::str_trunc(path_desc, 60)

            p_gsea_i <- gseaplot2(
                gsea_result,
                geneSetID = path_id,
                title     = path_desc,
                color     = ifelse(gsea_result@result$NES[
                    gsea_result@result$ID == path_id][1] > 0,
                    "#D73027", "#4575B4"),
                pvalue_table = TRUE
            )
            fn <- paste0("GSEA_enrichment_", i, "_", gsub("[^A-Za-z0-9]", "_", path_desc))
            fn <- substr(fn, 1, 80)
            ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_gsea_i,
                   width=12, height=7)
            ggsave(file.path(plot_dir, paste0(fn, ".png")), p_gsea_i,
                   width=12, height=7, dpi=300)
        }
    }, error=function(e) message("  Warning: enrichment plots failed: ", e$message))

    # 5. emapplot for GSEA
    message("  - emapplot")
    tryCatch({
        gsea_sim <- pairwise_termsim(gsea_result)
        p_emap <- emapplot(
            gsea_sim,
            showCategory = min(30, nrow(gsea_result)),
            color        = "NES",
            cex_label_category = 0.6
        ) +
            scale_color_gradient2(
                low="#4575B4", mid="grey90", high="#D73027",
                midpoint=0, name="NES"
            ) +
            labs(title="GSEA Enrichment Map",
                 subtitle="Color = NES; edge weight = gene overlap") +
            theme(plot.title = element_text(face="bold", size=14))
        save_plot(p_emap, "GSEA_emapplot", width=14, height=12)
    }, error=function(e) message("  Warning: GSEA emapplot failed: ", e$message))
}

# ── Session info ─────────────────────────────────────────────────────────────
sink(file.path(opt$outdir, "sessionInfo_gsea.txt"))
cat("=== GSEA - Session Information ===\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

message("\n=== GSEA analysis complete! ===")
