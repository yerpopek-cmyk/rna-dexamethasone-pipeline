#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis
# Project: RNA-seq Dexamethasone Effects (GSE52778)
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(EnhancedVolcano)
    library(RColorBrewer)
    library(dplyr)
    library(tibble)
    library(scales)
})

# ── Argument parsing ────────────────────────────────────────────────────────

option_list <- list(
    make_option("--counts",     type="character", help="featureCounts output file"),
    make_option("--metadata",   type="character", help="Sample metadata CSV"),
    make_option("--pvalue",     type="double",    default=0.05,  help="Adjusted p-value cutoff"),
    make_option("--lfc",        type="double",    default=1.0,   help="|log2FC| threshold"),
    make_option("--min_count",  type="integer",   default=10,    help="Minimum total count filter"),
    make_option("--lfc_shrink", type="logical",   default=TRUE,  help="Apply lfcShrink (apeglm)"),
    make_option("--outdir",     type="character", default=".",   help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

set.seed(42)

message("=== DESeq2 Analysis Pipeline ===")
message("Parameters:")
message("  Counts     : ", opt$counts)
message("  Metadata   : ", opt$metadata)
message("  p-value    : ", opt$pvalue)
message("  |LFC|      : ", opt$lfc)
message("  Min count  : ", opt$min_count)
message("  LFC shrink : ", opt$lfc_shrink)

# ── Load data ───────────────────────────────────────────────────────────────

message("\n[1/8] Loading count matrix and metadata...")

# Parse featureCounts output (skip comment lines starting with #)
counts_raw <- read.delim(opt$counts, comment.char="#", check.names=FALSE)

# Column structure: Geneid, Chr, Start, End, Strand, Length, sample1.bam, ...
gene_info  <- counts_raw[, 1:6]
count_mat  <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(count_mat) <- counts_raw$Geneid

# Clean sample names: remove path and .bam extension
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.bam$", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.sortedByCoord\\.out$", "", colnames(count_mat))

message("  Genes loaded   : ", nrow(count_mat))
message("  Samples loaded : ", ncol(count_mat))

# Load metadata
metadata <- read.csv(opt$metadata, stringsAsFactors=FALSE)
rownames(metadata) <- metadata$sample

# Align sample order
shared_samples <- intersect(colnames(count_mat), rownames(metadata))
if (length(shared_samples) == 0) {
    stop("ERROR: No matching samples between count matrix and metadata!\n",
         "Count matrix samples: ", paste(colnames(count_mat), collapse=", "), "\n",
         "Metadata samples: ",     paste(rownames(metadata), collapse=", "))
}
count_mat <- count_mat[, shared_samples]
metadata  <- metadata[shared_samples, , drop=FALSE]

message("  Matched samples: ", length(shared_samples), " (", 
        paste(shared_samples, collapse=", "), ")")

# ── Create DESeqDataSet ─────────────────────────────────────────────────────

message("\n[2/8] Creating DESeqDataSet...")

# Ensure condition is a factor with proper reference level
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref="control")

dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = metadata,
    design    = ~ condition
)

# ── Pre-filtering ────────────────────────────────────────────────────────────

message("\n[3/8] Filtering low-count genes...")
message("  Genes before filter: ", nrow(dds))

keep <- rowSums(counts(dds)) >= opt$min_count
dds  <- dds[keep, ]

message("  Genes after filter : ", nrow(dds))

# ── Run DESeq2 ───────────────────────────────────────────────────────────────

message("\n[4/8] Running DESeq2...")
dds <- DESeq(dds, parallel=FALSE)

# ── Extract results ──────────────────────────────────────────────────────────

message("\n[5/8] Extracting and shrinking results...")

contrast_levels <- levels(metadata$condition)
ref  <- contrast_levels[1]
test <- contrast_levels[2]
message("  Contrast: ", test, " vs ", ref)

# Raw results (for GSEA ranked list)
res_raw <- results(dds, 
                   contrast = c("condition", test, ref),
                   alpha    = opt$pvalue)

# Apply LFC shrinkage for publication plots
if (opt$lfc_shrink) {
    coef_name <- resultsNames(dds)[2]  # e.g. condition_treated_vs_control
    message("  Shrinkage coefficient: ", coef_name)
    if (requireNamespace("apeglm", quietly = TRUE)) {
        res <- lfcShrink(dds, coef=coef_name, type="apeglm")
    } else {
        message("  WARNING: 'apeglm' package not found. Falling back to 'normal' shrinkage method.")
        res <- lfcShrink(dds, coef=coef_name, type="normal")
    }
} else {
    res <- res_raw
}

# Convert to data frame
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj) %>%
    filter(!is.na(padj))

# Summary
n_up   <- sum(res_df$padj < opt$pvalue & res_df$log2FoldChange >  opt$lfc, na.rm=TRUE)
n_down <- sum(res_df$padj < opt$pvalue & res_df$log2FoldChange < -opt$lfc, na.rm=TRUE)
message("  Total DEGs (padj<", opt$pvalue, ", |LFC|>", opt$lfc, "): ", n_up + n_down)
message("    Up-regulated  : ", n_up)
message("    Down-regulated: ", n_down)

# Significant DEGs
res_sig <- res_df %>%
    filter(padj < opt$pvalue, abs(log2FoldChange) > opt$lfc)

# ── VST transformation ───────────────────────────────────────────────────────

message("\n[6/8] Computing VST for visualization...")
vsd <- vst(dds, blind=FALSE)

# ── Export tables ────────────────────────────────────────────────────────────

message("\n[7/8] Saving result tables...")

write.csv(res_df, 
          file.path(opt$outdir, "DEG_results.csv"), 
          row.names=FALSE, quote=FALSE)

write.csv(res_sig, 
          file.path(opt$outdir, "DEG_significant.csv"), 
          row.names=FALSE, quote=FALSE)

write.csv(as.data.frame(counts(dds, normalized=TRUE)) %>% rownames_to_column("gene_id"),
          file.path(opt$outdir, "normalized_counts.csv"),
          row.names=FALSE, quote=FALSE)

write.csv(as.data.frame(assay(vsd)) %>% rownames_to_column("gene_id"),
          file.path(opt$outdir, "vst_counts.csv"),
          row.names=FALSE, quote=FALSE)

saveRDS(dds, file.path(opt$outdir, "rds", "dds.rds"))
saveRDS(res, file.path(opt$outdir, "rds", "res.rds"))
saveRDS(vsd, file.path(opt$outdir, "rds", "vsd.rds"))

# ── Plots ────────────────────────────────────────────────────────────────────

message("\n[8/8] Generating plots...")

plot_dir <- file.path(opt$outdir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

save_plot <- function(plot_obj, filename, width=10, height=8) {
    base <- file.path(plot_dir, filename)
    ggsave(paste0(base, ".pdf"), plot_obj, width=width, height=height, dpi=300)
    ggsave(paste0(base, ".png"), plot_obj, width=width, height=height, dpi=300)
}

# ── PCA plot ─────────────────────────────────────────────────────────────────
message("  - PCA plot")
pca_data <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_labels <- as.data.frame(colData(dds))
pca_data$cell <- pca_labels$cell_line[match(rownames(pca_labels), pca_data$name)]

p_pca <- ggplot(pca_data, aes(PC1, PC2, color=condition, label=name)) +
    geom_point(size=5, alpha=0.9) +
    ggrepel::geom_text_repel(size=3.5, show.legend=FALSE) +
    scale_color_manual(values=c("control"="#4575B4", "treated"="#D73027")) +
    labs(
        title    = "Principal Component Analysis",
        subtitle = paste("Variance stabilized counts (n =", nrow(dds), "genes)"),
        x        = paste0("PC1: ", pct_var[1], "% variance"),
        y        = paste0("PC2: ", pct_var[2], "% variance"),
        color    = "Condition"
    ) +
    theme_bw(base_size=14) +
    theme(
        plot.title   = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        panel.grid.minor = element_blank()
    )

save_plot(p_pca, "PCA_plot")

# ── MA plot ──────────────────────────────────────────────────────────────────
message("  - MA plot")
res_ma <- as.data.frame(res_raw) %>%
    rownames_to_column("gene_id") %>%
    mutate(
        significant = !is.na(padj) & padj < opt$pvalue & abs(log2FoldChange) > opt$lfc,
        label_gene  = ifelse(significant & abs(log2FoldChange) > 3, gene_id, NA)
    )

p_ma <- ggplot(res_ma, aes(x=log10(baseMean + 1), y=log2FoldChange)) +
    geom_point(aes(color=significant), alpha=0.4, size=0.8) +
    geom_hline(yintercept=c(-opt$lfc, opt$lfc), linetype="dashed", color="grey40") +
    geom_hline(yintercept=0, color="black") +
    scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#E41A1C"),
                       labels=c("Not significant", "Significant DEG")) +
    labs(
        title    = "MA Plot: Differential Expression",
        subtitle = paste0("Dexamethasone treated vs Control (padj<",
                          opt$pvalue, ", |LFC|>", opt$lfc, ")"),
        x        = "log10(Mean normalized counts + 1)",
        y        = "log2 Fold Change",
        color    = NULL
    ) +
    theme_bw(base_size=14) +
    theme(
        plot.title       = element_text(face="bold"),
        legend.position  = "top",
        panel.grid.minor = element_blank()
    )

save_plot(p_ma, "MA_plot")

# ── Volcano plot ─────────────────────────────────────────────────────────────
message("  - Volcano plot")

# Label top 30 most significant genes
top_labels <- res_df %>%
    filter(!is.na(padj)) %>%
    slice_min(order_by=padj, n=30) %>%
    pull(gene_id)

keyvals <- ifelse(
    res_df$log2FoldChange >  opt$lfc & res_df$padj < opt$pvalue, "#D73027",
    ifelse(
        res_df$log2FoldChange < -opt$lfc & res_df$padj < opt$pvalue, "#4575B4",
        "grey60"
    )
)
names(keyvals) <- ifelse(
    res_df$log2FoldChange >  opt$lfc & res_df$padj < opt$pvalue, "Up-regulated",
    ifelse(
        res_df$log2FoldChange < -opt$lfc & res_df$padj < opt$pvalue, "Down-regulated",
        "Not significant"
    )
)

ev_labels <- ifelse(res_df$gene_id %in% top_labels, res_df$gene_id, NA)

p_volcano <- EnhancedVolcano(
    res_df,
    lab           = ev_labels,
    x             = "log2FoldChange",
    y             = "padj",
    title         = "Volcano Plot: Dexamethasone vs Control",
    subtitle      = paste("Up:", n_up, "  |  Down:", n_down,
                          "  |  p.adj <", opt$pvalue, ", |LFC| >", opt$lfc),
    pCutoff       = opt$pvalue,
    FCcutoff      = opt$lfc,
    pointSize     = 1.5,
    labSize       = 3.0,
    colCustom     = keyvals,
    legendPosition= "right",
    legendLabSize = 10,
    axisLabSize   = 12,
    titleLabSize  = 14,
    col           = c("grey60", "#4575B4", "#D73027", "darkred"),
    gridlines.major = TRUE,
    gridlines.minor = FALSE
)

save_plot(p_volcano, "Volcano_plot", width=12, height=10)

# ── Dispersion plot ──────────────────────────────────────────────────────────
message("  - Dispersion estimates plot")
png(file.path(plot_dir, "Dispersion_plot.png"), width=2400, height=1800, res=300)
plotDispEsts(dds, main="DESeq2 Dispersion Estimates")
dev.off()

# ── Heatmap top DEGs ─────────────────────────────────────────────────────────
message("  - Heatmap top 50 DEGs")
n_top  <- min(50, nrow(res_sig))
top_de <- head(res_sig$gene_id, n_top)
mat    <- assay(vsd)[top_de, ]

# Z-score by row for better visualization
mat_scaled <- t(scale(t(mat)))

annotation_col <- as.data.frame(colData(dds))[, "condition", drop=FALSE]
ann_colors <- list(
    condition = c(control="#4575B4", treated="#D73027")
)

pheatmap_plot <- pheatmap(
    mat_scaled,
    annotation_col   = annotation_col,
    annotation_colors= ann_colors,
    show_rownames    = n_top <= 30,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks           = seq(-3, 3, length.out=101),
    main             = paste0("Top ", n_top, " DEGs — Z-scored VST counts"),
    fontsize         = 10,
    fontsize_row     = 7,
    border_color     = NA,
    silent           = TRUE
)

ggsave(file.path(plot_dir, "Heatmap_top_DEGs.pdf"),
       pheatmap_plot$gtable, width=10, height=14)
ggsave(file.path(plot_dir, "Heatmap_top_DEGs.png"),
       pheatmap_plot$gtable, width=10, height=14, dpi=300)

# ── Sample correlation heatmap ────────────────────────────────────────────────
message("  - Sample correlation heatmap")
vst_mat <- assay(vsd)
sampleDists <- dist(t(vst_mat))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

corr_heatmap <- pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    annotation_col  = annotation_col,
    annotation_colors= ann_colors,
    color           = colorRampPalette(c("#4575B4", "white"))(100),
    main            = "Sample-to-sample Euclidean distances (VST)",
    fontsize        = 11,
    border_color    = NA,
    silent          = TRUE
)

ggsave(file.path(plot_dir, "Sample_correlation_heatmap.pdf"),
       corr_heatmap$gtable, width=8, height=7)
ggsave(file.path(plot_dir, "Sample_correlation_heatmap.png"),
       corr_heatmap$gtable, width=8, height=7, dpi=300)

# ── Count plot for top DEGs ──────────────────────────────────────────────────
message("  - Normalized counts for top 9 DEGs")
top9_genes <- head(res_sig$gene_id, 9)

count_df <- lapply(top9_genes, function(gene) {
    cnt <- plotCounts(dds, gene=gene, intgroup="condition", returnData=TRUE)
    cnt$gene <- gene
    cnt
}) %>% bind_rows()

p_counts <- ggplot(count_df, aes(x=condition, y=count, color=condition)) +
    geom_jitter(width=0.15, size=2.5, alpha=0.8) +
    geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.4) +
    scale_y_log10(labels=scales::comma) +
    scale_color_manual(values=c("control"="#4575B4", "treated"="#D73027")) +
    facet_wrap(~ gene, scales="free_y", ncol=3) +
    labs(
        title = "Normalized Counts for Top Differentially Expressed Genes",
        x     = "Condition",
        y     = "Normalized count (log10 scale)",
        color = "Condition"
    ) +
    theme_bw(base_size=11) +
    theme(
        plot.title   = element_text(face="bold"),
        strip.text   = element_text(face="bold.italic"),
        legend.position = "none"
    )

save_plot(p_counts, "Top9_DEG_counts", width=12, height=10)

# ── Session info ─────────────────────────────────────────────────────────────
sink(file.path(opt$outdir, "sessionInfo_deseq2.txt"))
cat("=== DESeq2 Analysis - Session Information ===\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

message("\n=== DESeq2 analysis complete! ===")
message("Results saved to: ", opt$outdir)
message("  DEGs (total)  : ", n_up + n_down, " (", n_up, " up, ", n_down, " down)")
