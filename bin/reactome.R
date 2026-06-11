#!/usr/bin/env Rscript
# =============================================================================
# Reactome Pathway Enrichment Analysis (ORA)
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
    library(ggupset)
    library(cowplot)
})

# ── Argument parsing ────────────────────────────────────────────────────────

option_list <- list(
    make_option("--deg",      type="character", help="DEG significant CSV (from deseq2.R)"),
    make_option("--organism", type="character", default="human", help="Organism for ReactomePA"),
    make_option("--pvalue",   type="double",    default=0.05,    help="p-value cutoff"),
    make_option("--lfc",      type="double",    default=1.0,     help="|log2FC| threshold"),
    make_option("--outdir",   type="character", default=".",     help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

set.seed(42)

message("=== Reactome Pathway Enrichment (ORA) ===")
message("Parameters:")
message("  DEG input   : ", opt$deg)
message("  Organism    : ", opt$organism)
message("  p-value     : ", opt$pvalue)

plot_dir <- file.path(opt$outdir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

save_plot <- function(plot_obj, filename, width=12, height=10) {
    base <- file.path(plot_dir, filename)
    tryCatch({
        ggsave(paste0(base, ".pdf"), plot_obj, width=width, height=height)
        ggsave(paste0(base, ".png"), plot_obj, width=width, height=height, dpi=300)
    }, error=function(e) message("  Warning: could not save ", filename, ": ", e$message))
}

# ── Load DEGs ───────────────────────────────────────────────────────────────

message("\n[1/5] Loading significant DEGs...")
deg <- read.csv(opt$deg)
if (nrow(deg) > 0 && "gene_id" %in% colnames(deg)) {
    deg$gene_id <- sub("\\..*$", "", deg$gene_id)
}
message("  Loaded DEGs: ", nrow(deg))

gene_ids <- deg$gene_id

# ── ID Conversion ────────────────────────────────────────────────────────────

message("\n[2/5] Converting Ensembl IDs to Entrez IDs...")

# Determine ID type (Ensembl vs gene symbol)
is_ensembl <- grepl("^ENSG", gene_ids[1])

if (is_ensembl) {
    from_type <- "ENSEMBL"
} else {
    from_type <- "SYMBOL"
}

entrez_df <- bitr(
    gene_ids,
    fromType = from_type,
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db,
    drop     = TRUE
)

# Also get gene symbols for readable output
symbol_df <- bitr(
    entrez_df$ENTREZID,
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = org.Hs.eg.db,
    drop     = TRUE
)

entrez_ids <- unique(entrez_df$ENTREZID)
message("  Input genes      : ", length(gene_ids))
message("  Successfully mapped: ", length(entrez_ids))
message("  Mapping rate     : ", round(length(entrez_ids)/length(gene_ids)*100, 1), "%")

if (length(entrez_ids) < 5) {
    stop("ERROR: Too few genes mapped to Entrez IDs (", length(entrez_ids), "). ",
         "Check that the gene_id column contains valid Ensembl or HGNC symbols.")
}

# ── Separate up/down regulated genes ─────────────────────────────────────────

up_genes   <- deg %>% filter(log2FoldChange >  0) %>% pull(gene_id)
down_genes <- deg %>% filter(log2FoldChange <  0) %>% pull(gene_id)

up_entrez   <- bitr(up_genes,   from_type, "ENTREZID", org.Hs.eg.db, drop=TRUE)$ENTREZID
down_entrez <- bitr(down_genes, from_type, "ENTREZID", org.Hs.eg.db, drop=TRUE)$ENTREZID

# ── Reactome ORA ─────────────────────────────────────────────────────────────

message("\n[3/5] Running Reactome Over-Representation Analysis...")

run_enrichPathway <- function(gene_list, label) {
    if (length(gene_list) < 3) {
        message("  Skipping ", label, ": fewer than 3 genes")
        return(NULL)
    }
    result <- tryCatch(
        enrichPathway(
            gene         = gene_list,
            organism     = opt$organism,
            pvalueCutoff = opt$pvalue,
            pAdjustMethod= "BH",
            qvalueCutoff = 0.2,
            readable     = TRUE,
            minGSSize    = 10,
            maxGSSize    = 500
        ),
        error = function(e) {
            message("  Warning: enrichPathway failed for ", label, ": ", e$message)
            return(NULL)
        }
    )
    if (!is.null(result) && nrow(result) > 0) {
        message("  ", label, ": found ", nrow(result), " enriched pathways")
    } else {
        message("  ", label, ": no significant pathways found")
    }
    return(result)
}

enrich_all  <- run_enrichPathway(entrez_ids, "All DEGs")
enrich_up   <- run_enrichPathway(up_entrez,   "Up-regulated")
enrich_down <- run_enrichPathway(down_entrez,  "Down-regulated")

# ── Export results ────────────────────────────────────────────────────────────

message("\n[4/5] Saving results...")

save_result <- function(enrich_obj, filename) {
    if (!is.null(enrich_obj) && nrow(enrich_obj) > 0) {
        df <- as.data.frame(enrich_obj)
        write.csv(df, file.path(opt$outdir, filename), row.names=FALSE, quote=FALSE)
        message("  Saved: ", filename, " (", nrow(df), " pathways)")
        return(df)
    }
    message("  No results for: ", filename)
    return(NULL)
}

df_all  <- save_result(enrich_all,  "reactome_ORA_results.csv")
df_up   <- save_result(enrich_up,   "reactome_ORA_upregulated.csv")
df_down <- save_result(enrich_down, "reactome_ORA_downregulated.csv")

# ── Plots ─────────────────────────────────────────────────────────────────────

message("\n[5/5] Generating pathway plots...")

if (!is.null(enrich_all) && nrow(enrich_all) > 0) {

    n_show <- min(20, nrow(enrich_all))

    # 1. Dotplot
    message("  - Dotplot")
    p_dot <- dotplot(enrich_all, showCategory=n_show) +
        labs(title="Reactome ORA: Top Enriched Pathways",
             subtitle=paste("All DEGs | padj <", opt$pvalue)) +
        theme(
            plot.title   = element_text(face="bold", size=14),
            axis.text.y  = element_text(size=9),
            legend.position = "right"
        )
    save_plot(p_dot, "ORA_dotplot", width=14, height=10)

    # 2. Barplot
    message("  - Barplot")
    p_bar <- barplot(enrich_all, showCategory=n_show) +
        labs(title="Reactome ORA: Pathway Gene Count") +
        theme(
            plot.title  = element_text(face="bold", size=14),
            axis.text.y = element_text(size=9)
        )
    save_plot(p_bar, "ORA_barplot", width=14, height=10)

    # 3. cnetplot (gene-concept network)
    message("  - cnetplot")
    tryCatch({
        fold_changes <- setNames(deg$log2FoldChange, deg$gene_id)
        fc_entrez    <- setNames(
            fold_changes[entrez_df[[from_type]]],
            entrez_df$ENTREZID
        )
        fc_entrez <- fc_entrez[!is.na(fc_entrez)]

        p_cnet <- cnetplot(
            enrich_all,
            foldChange   = fc_entrez,
            categorySize = "pvalue",
            showCategory = 6,
            circular     = FALSE,
            colorEdge    = TRUE,
            node_label   = "all"
        ) +
            scale_color_gradient2(
                low      = "#4575B4",
                mid      = "grey90",
                high     = "#D73027",
                midpoint = 0,
                name     = "log2FC"
            ) +
            labs(title="Gene-Pathway Network",
                 subtitle="Node size = gene ratio; edge = gene-pathway link") +
            theme(plot.title = element_text(face="bold", size=14))

        save_plot(p_cnet, "ORA_cnetplot", width=14, height=12)
    }, error=function(e) message("  Warning: cnetplot failed: ", e$message))

    # 4. emapplot (enrichment map)
    message("  - emapplot")
    tryCatch({
        enrich_sim <- pairwise_termsim(enrich_all)
        p_emap <- emapplot(
            enrich_sim,
            showCategory = min(30, nrow(enrich_all)),
            color        = "p.adjust",
            cex_label_category = 0.6
        ) +
            labs(title="Enrichment Map: Pathway Similarity",
                 subtitle="Edge weight = gene overlap between pathways") +
            theme(plot.title = element_text(face="bold", size=14))

        save_plot(p_emap, "ORA_emapplot", width=14, height=12)
    }, error=function(e) message("  Warning: emapplot failed: ", e$message))

    # 5. Upsetplot
    message("  - Upsetplot")
    tryCatch({
        p_upset <- upsetplot(enrich_all) +
            labs(title="Gene Membership Across Pathways")
        save_plot(p_upset, "ORA_upsetplot", width=12, height=8)
    }, error=function(e) message("  Warning: upsetplot failed: ", e$message))

    # 6. Up vs Down comparison plot
    if (!is.null(enrich_up) && !is.null(enrich_down) &&
        nrow(enrich_up) > 0 && nrow(enrich_down) > 0) {
        message("  - Up vs Down comparison dotplot")
        tryCatch({
            compare_list <- list(Up=up_entrez, Down=down_entrez)
            compare_result <- compareCluster(
                geneClusters = compare_list,
                fun          = "enrichPathway",
                organism     = opt$organism,
                pvalueCutoff = opt$pvalue,
                readable     = TRUE
            )
            p_compare <- dotplot(compare_result, showCategory=15) +
                scale_color_gradient(low="#D73027", high="#4575B4") +
                labs(title="Pathway Enrichment: Up vs Down-regulated Genes") +
                theme(
                    plot.title  = element_text(face="bold", size=14),
                    axis.text.y = element_text(size=8)
                )
            save_plot(p_compare, "ORA_up_vs_down_dotplot", width=14, height=12)
        }, error=function(e) message("  Warning: compareCluster plot failed: ", e$message))
    }
}

# ── Session info ─────────────────────────────────────────────────────────────
sink(file.path(opt$outdir, "sessionInfo_reactome.txt"))
cat("=== Reactome ORA - Session Information ===\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

message("\n=== Reactome ORA complete! ===")
