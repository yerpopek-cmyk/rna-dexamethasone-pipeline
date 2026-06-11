<div align="center">

# 🧬 rna-dexamethasone

### Production-grade RNA-seq Pipeline
**Dexamethasone Effects on Human Airway Smooth Muscle Cells**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![R 4.5+](https://img.shields.io/badge/R-%E2%89%A54.5-276DC3?logo=R)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

---

## Table of Contents

1. [What this project is about](#what-this-project-is-about)
2. [The biological question](#the-biological-question)
3. [Dataset — GSE52778](#dataset--gse52778)
4. [Pipeline architecture](#pipeline-architecture)
5. [Step-by-step walkthrough](#step-by-step-walkthrough)
   - [Step 1 — Quality control (FastQC)](#step-1--quality-control-fastqc)
   - [Step 2 — Adapter trimming (fastp)](#step-2--adapter-trimming-fastp)
   - [Step 3 — Genome alignment (STAR)](#step-3--genome-alignment-star)
   - [Step 4 — BAM quality control (samtools)](#step-4--bam-quality-control-samtools)
   - [Step 5 — Read quantification (featureCounts)](#step-5--read-quantification-featurecounts)
   - [Step 6 — Differential expression (DESeq2)](#step-6--differential-expression-deseq2)
   - [Step 7 — Pathway analysis — ORA (ReactomePA)](#step-7--pathway-analysis--ora-reactomepa)
   - [Step 8 — Pathway analysis — GSEA (ReactomePA)](#step-8--pathway-analysis--gsea-reactomepa)
   - [Step 9 — HTML report (RMarkdown)](#step-9--html-report-rmarkdown)
   - [Step 10 — Quality summary (MultiQC)](#step-10--quality-summary-multiqc)
6. [Understanding every output file](#understanding-every-output-file)
7. [Repository structure](#repository-structure)
8. [Requirements](#requirements)
9. [Installation](#installation)
10. [Running the pipeline](#running-the-pipeline)
11. [Parameters reference](#parameters-reference)
12. [Samplesheet format](#samplesheet-format)
13. [Expected results](#expected-results-gse52778)
14. [Troubleshooting](#troubleshooting)
15. [Citations](#citations)
16. [License](#license)

---

## What this project is about

This repository contains a **complete, end-to-end RNA-sequencing analysis pipeline** built with [Nextflow DSL2](https://www.nextflow.io/) and R/Bioconductor. Starting from raw FASTQ files downloaded from the public [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/), the pipeline produces:

- A genome-aligned BAM file per sample with full quality metrics
- A gene-level read count matrix (genes × samples)
- A table of differentially expressed genes (DEGs) with statistics
- Publication-ready plots: volcano plot, PCA, heatmaps, MA plot
- Reactome pathway enrichment results (ORA and GSEA)
- A single self-contained HTML report summarising everything
- A MultiQC dashboard covering all quality steps

Everything runs inside a **Docker container** that packages every tool and R library needed — you do not need to install anything except Nextflow and Docker.

---

## The biological question

> **How does dexamethasone change gene expression in human airway smooth muscle cells?**

Dexamethasone (Dex) is a potent synthetic glucocorticoid used to treat asthma, COPD, and inflammatory diseases. It works primarily by binding the glucocorticoid receptor (GR), which then enters the nucleus and alters gene transcription on a genome-wide scale.

The key outcomes we measure:

| Question | Method |
|----------|--------|
| Which genes are up/down regulated? | DESeq2 differential expression |
| Are the changes statistically sound? | Negative binomial modelling, BH-adjusted p-values |
| Which biological pathways are affected? | Reactome ORA (gene-list enrichment) |
| Does the full transcriptome confirm pathway activation/suppression? | Reactome GSEA (ranked list enrichment) |

---

## Dataset — GSE52778

| Property | Value |
|----------|-------|
| GEO Accession | [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) |
| Publication | Himes et al., PLOS ONE, 2014 |
| Organism | *Homo sapiens* |
| Cell type | Primary human airway smooth muscle (HASM) cells |
| Treatment | 1 µM dexamethasone for 18 hours |
| Design | 4 cell lines (donors), paired control + treated = **8 samples** |
| Sequencing | Paired-end 75 bp, Illumina HiSeq 2000 |
| Reference | GRCh38 (GENCODE v46) |

### The 8 samples

| SRR ID | Condition | Cell Line | Batch | Sex |
|--------|-----------|-----------|-------|-----|
| SRR1039508 | Control | N61311 | 1 | Male |
| SRR1039509 | Treated | N61311 | 1 | Male |
| SRR1039512 | Control | N052611 | 1 | Male |
| SRR1039513 | Treated | N052611 | 1 | Male |
| SRR1039516 | Control | N080611 | 2 | Female |
| SRR1039517 | Treated | N080611 | 2 | Female |
| SRR1039520 | Control | N061011 | 2 | Female |
| SRR1039521 | Treated | N061011 | 2 | Female |

The paired design (same donor, treated vs untreated) is a major strength — it eliminates inter-donor variability, making DESeq2 highly powered even with only 4 biological replicates per condition.

---

## Pipeline architecture

```
FASTQ files (8 samples, paired-end)
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│  Subworkflow: QC_TRIM                                       │
│  ┌─────────────┐    ┌────────────────────────────────────┐  │
│  │   FASTQC    │    │              FASTP                 │  │
│  │ (raw reads) │    │  (trim adapters, low-quality ends) │  │
│  └─────────────┘    └────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
         │ trimmed FASTQ (8 samples)
         ▼
┌─────────────────────────────────────────────────────────────┐
│  Subworkflow: ALIGN_QUANTIFY                                │
│  ┌─────────────────────────────────────────────────────┐    │
│  │  STAR_ALIGN  (map to GRCh38, produce sorted BAMs)   │    │
│  └─────────────────────────────────────────────────────┘    │
│  ┌─────────────────────────────────────────────────────┐    │
│  │  BAM_QC  (flagstat, stats, idxstats per BAM)        │    │
│  └─────────────────────────────────────────────────────┘    │
│  ┌─────────────────────────────────────────────────────┐    │
│  │  FEATURECOUNTS  (gene-level count matrix, all 8)    │    │
│  └─────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────┘
         │ gene_counts.txt  +  metadata.csv
         ▼
┌─────────────────────────────────────────────────────────────┐
│  Subworkflow: DIFF_EXPRESSION                               │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  DESEQ2  → DEGs, PCA, volcano, heatmap, MA plot      │   │
│  └──────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  REACTOME_ORA  → pathway enrichment (DEG list)       │   │
│  └──────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  REACTOME_GSEA → pathway enrichment (ranked list)    │   │
│  └──────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  RENDER_REPORT → RNAseq_analysis_report.html         │   │
│  └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
         │ (in parallel with everything above)
         ▼
┌─────────────────────────────────────────────────────────────┐
│  MODULE: MULTIQC                                            │
│  (aggregates FastQC + fastp + STAR + featureCounts + BAM    │
│   QC logs into a single interactive HTML dashboard)         │
└─────────────────────────────────────────────────────────────┘
```

---

## Step-by-step walkthrough

### Step 1 — Quality control (FastQC)

**Module:** [`modules/fastqc.nf`](modules/fastqc.nf)  
**Tool:** FastQC v0.12.1  
**What it does:**

FastQC inspects each raw FASTQ file and generates a quality report. It checks ~11 metrics including:

- **Per-base sequence quality** — Phred scores per cycle. Good reads have median Q ≥ 28 across all positions.
- **Per-sequence quality scores** — Distribution of mean quality per read. Should peak above Q30.
- **Per-base sequence content** — Base composition per cycle. Should be flat (random) after the first few bases.
- **GC content** — Should follow a smooth normal distribution centred on the organism's expected GC (human ~50%).
- **Sequence duplication levels** — Very high duplication can indicate PCR artefacts.
- **Adapter content** — Illumina TruSeq adapters should be absent or minimal in good data; fastp will trim them.

**Outputs** (in `results/fastqc/`):
- `*_fastqc.html` — Interactive report (open in browser)
- `*_fastqc.zip` — Raw data for MultiQC

> **Why this matters:** If a sample fails FastQC (e.g., GC spikes, heavy adapter contamination), you know before wasting compute time on alignment.

---

### Step 2 — Adapter trimming (fastp)

**Module:** [`modules/fastp.nf`](modules/fastp.nf)  
**Tool:** fastp v0.23.4  
**What it does:**

fastp removes Illumina adapter sequences and low-quality bases from each read end. Default settings:
- Adapter auto-detection (detects TruSeq automatically)
- Quality window trimming: trim when sliding-window Q < 20 (`--qualified_quality_phred 20`)
- Minimum read length: 30 bp after trimming (`--length_required 30`)
- Both R1 and R2 trimmed together, preserving read pairing

**Outputs** (in `results/fastp/`):
- `*.fastp.html` — Per-sample trimming report
- `*.fastp.json` — Machine-readable statistics
- `*.fastp.log` — Parsed by MultiQC

> **Why this matters:** Adapter sequences in reads cause them to fail alignment. Low-quality tails reduce mapping accuracy. Trimming is a mandatory pre-processing step.

---

### Step 3 — Genome alignment (STAR)

**Module:** [`modules/star_align.nf`](modules/star_align.nf)  
**Tool:** STAR v2.7.11b  
**Reference:** GRCh38 primary assembly + GENCODE v46 annotation  
**What it does:**

STAR is a splice-aware aligner. "Splice-aware" is critical for RNA-seq: reads that span exon-exon junctions would fail to align with DNA-centric aligners. STAR maps each read pair against the pre-built genome index and:

1. Finds the longest prefix of a read that matches uniquely in the genome
2. Identifies the splice junction and extends the alignment to the next exon
3. Outputs a coordinate-sorted BAM file per sample

Key flags used:
- `--outSAMtype BAM SortedByCoordinate` — produces a position-sorted BAM (required by featureCounts and most downstream tools)
- `--outSAMattributes NH HI AS NM MD` — SAM attributes for multi-mapping (NH), alignment score (AS), and mismatch info (NM, MD)
- `--outFilterIntronMotifs RemoveNoncanonical` — discards reads aligned via non-canonical splice sites
- `--limitBAMsortRAM 1500000000` — caps BAM sorting RAM at 1.5 GB (safe for most laptops)
- `--runThreadN` — multi-threaded (4 threads by default)

After alignment, `samtools index` creates a `.bai` index file so any tool can perform random access into the BAM.

**Outputs** (in `results/star/<sample_id>/`):
- `bam/*.sortedByCoord.out.bam` — Coordinate-sorted aligned reads
- `bam/*.sortedByCoord.out.bam.bai` — BAM index
- `log/*Log.final.out` — Summary of mapping statistics (% mapped, multi-mappers, etc.)
- `spliceJunctions/*.SJ.out.tab` — All splice junctions detected in this sample

> **Why this matters:** The fraction of reads that map uniquely to the genome directly determines analysis power. Expected: > 90% uniquely mapped for human samples.

---

### Step 4 — BAM quality control (samtools)

**Module:** [`modules/bam_qc.nf`](modules/bam_qc.nf)  
**Tool:** samtools (bundled in Docker image)  
**What it does:**

Runs three complementary QC commands on each BAM file:

| Command | Output | What it shows |
|---------|--------|---------------|
| `samtools flagstat` | `*.flagstat` | Total reads, mapped, paired, properly paired, duplicates |
| `samtools stats` | `*.stats` | Insert size, base quality, alignment error rates, coverage stats |
| `samtools idxstats` | `*.idxstats` | Read counts per chromosome (checks for rDNA/mitochondrial contamination) |

**Outputs** (in `results/bam_qc/`):
- `*.flagstat` — One-line summary counts per sample
- `*.stats` — Detailed statistics (MultiQC parses ~20 metrics from this)
- `*.idxstats` — Per-chromosome read distribution

> **Why this matters:** A properly paired rate < 80% or an unusually high mitochondrial fraction (MT chromosome in idxstats) signals a problem with sample quality.

---

### Step 5 — Read quantification (featureCounts)

**Module:** [`modules/featurecounts.nf`](modules/featurecounts.nf)  
**Tool:** Subread featureCounts v2.0.6  
**What it does:**

featureCounts counts how many reads (read pairs) overlap each annotated gene in the GTF. It runs on **all 8 BAM files simultaneously**, producing a single count matrix.

Key flags:
- `-t exon` — count reads overlapping exon features (not introns)
- `-g gene_id` — summarise to gene level (not transcript level)
- `-B` — require both reads in a pair to be aligned (for paired-end data)
- `-C` — do not count read pairs whose two reads map to different chromosomes (chimeric)
- `--minOverlap 10` — a read must overlap an exon by at least 10 bp to be counted
- `-p --countReadPairs` — count read pairs (fragments), not individual reads

**Outputs** (in `results/featurecounts/`):
- `gene_counts.txt` — The count matrix: rows = ENSEMBL gene IDs, columns = samples. This is the primary input to DESeq2.
- `gene_counts.txt.summary` — Assignment statistics: how many reads were assigned, ambiguous, unmapped, etc.

> **Why this matters:** The count matrix is the bridge between sequencing and statistics. Genes with > 10 total counts across all samples are kept for DE analysis.

---

### Step 6 — Differential expression (DESeq2)

**Script:** [`bin/deseq2.R`](bin/deseq2.R)  
**Tool:** DESeq2 v1.40+, apeglm  
**What it does (8 steps internally):**

#### 6.1 — Load data
Reads `gene_counts.txt` from featureCounts. Strips `.sortedByCoord.out.bam` suffixes from column names to match sample names in `metadata.csv`.

#### 6.2 — Create DESeqDataSet
Builds the DESeq2 data object with design `~ condition` (comparing treated vs. control). The condition factor is releveled so "control" is the reference — meaning all fold changes are reported as **treated ÷ control**.

#### 6.3 — Pre-filtering
Removes genes where the total count across all 8 samples is < 10. This eliminates ~35,000 uninformative genes (mostly lncRNAs and pseudogenes with no detected expression), reducing the multiple testing burden and speeding up computation.

#### 6.4 — Run DESeq2
Fits a negative binomial generalised linear model per gene. The negative binomial distribution is used because RNA-seq counts are over-dispersed (variance > mean). DESeq2 estimates dispersion with shrinkage towards a common trend — this is what makes it robust with small sample sizes.

#### 6.5 — LFC shrinkage (apeglm)
Applies **apeglm** shrinkage to log2 fold changes. Without shrinkage, genes with very low counts produce artificially extreme fold changes (e.g., 0 counts in control, 1 count in treated = log2FC = ∞). apeglm shrinks these estimates towards zero while preserving large, confident fold changes. This is the recommended method for volcano plots and ranking.

#### 6.6 — VST transformation
Computes the Variance Stabilising Transformation (VST). VST converts raw counts into a scale where the variance is approximately independent of the mean — making distance-based methods (PCA, clustering, heatmaps) meaningful. It is analogous to log2(CPM) but statistically principled.

#### 6.7 — Export tables
- `DEG_results.csv` — All tested genes with: `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`
- `DEG_significant.csv` — Subset: padj < 0.05 AND |log2FC| > 1 (i.e., ≥ 2-fold change)
- `normalized_counts.csv` — DESeq2 size-factor normalised counts (corrects for library size)
- `vst_counts.csv` — VST-transformed matrix for downstream visualisation

#### 6.8 — Generate plots

| Plot | What you see | What it means |
|------|-------------|---------------|
| `PCA_plot` | 8 samples as dots in 2D; colour = condition | PC1 should separate control from treated if there is a strong treatment effect |
| `MA_plot` | log2FC (y) vs. mean expression (x); red = significant | Low-count genes have scattered fold changes; significant genes (red) should form a cloud |
| `Volcano_plot` | -log10(padj) (y) vs. log2FC (x); red = up, blue = down | Top 30 genes labelled; shows the number and magnitude of DEGs |
| `Heatmap_top_DEGs` | Z-scored VST counts for top 50 DEGs; clustered | Clear block pattern confirms control vs. treated separation |
| `Sample_correlation_heatmap` | Euclidean distances between samples in VST space | Samples from the same cell line should cluster together |
| `Top9_DEG_counts` | Jitter + boxplot of normalised counts for top 9 genes | Shows raw effect size and within-group variability per gene |
| `Dispersion_plot` | Gene-wise dispersion estimates vs. mean | Should show the characteristic "V" pattern converging to the fitted trend |

---

### Step 7 — Pathway analysis — ORA (ReactomePA)

**Script:** [`bin/reactome.R`](bin/reactome.R)  
**Tool:** ReactomePA, clusterProfiler, org.Hs.eg.db  
**Method:** Over-Representation Analysis (ORA)

#### What is ORA?

ORA asks: *"Given my list of significant DEGs, are certain biological pathways more represented than expected by chance?"*

It uses a hypergeometric test (Fisher's exact test equivalent) for each Reactome pathway:
- **Universe**: all tested genes (~22,000 after filtering)
- **Gene list**: significant DEGs from DESeq2 (padj < 0.05, |LFC| > 1)
- **Pathway gene set**: genes annotated to each Reactome pathway

A pathway is significant if its BH-adjusted p-value < 0.05.

#### Steps internally:

1. **Load significant DEGs** — reads `DEG_significant.csv`
2. **ID conversion** — Ensembl gene IDs (ENSG...) → Entrez IDs using `org.Hs.eg.db`. This is required because ReactomePA uses Entrez IDs internally.
3. **Separate up/down genes** — runs ORA three times: all DEGs, up-regulated only, down-regulated only
4. **enrichPathway()** — the core test; `minGSSize=10`, `maxGSSize=500`, `qvalueCutoff=0.2`
5. **Save results** + **generate plots**

#### ORA output plots

| Plot | What it shows |
|------|--------------|
| `ORA_dotplot` | Top 20 pathways; dot size = gene ratio; colour = p.adjust |
| `ORA_barplot` | Top 20 pathways; bar length = gene count in pathway |
| `ORA_cnetplot` | Network: pathways (large nodes) connected to their member genes (small nodes); gene colour = log2FC |
| `ORA_emapplot` | Pathway similarity network; edges = shared genes; helps identify pathway clusters |
| `ORA_upsetplot` | Shows which genes are shared across multiple pathways (intersection matrix) |
| `ORA_up_vs_down_dotplot` | Side-by-side comparison of pathways enriched in up vs. down-regulated genes |

---

### Step 8 — Pathway analysis — GSEA (ReactomePA)

**Script:** [`bin/gsea.R`](bin/gsea.R)  
**Tool:** ReactomePA::gsePathway, fgsea  
**Method:** Gene Set Enrichment Analysis (GSEA)

#### What is GSEA and how is it different from ORA?

ORA uses only a hard-filtered gene list (requires an arbitrary threshold). GSEA uses **all tested genes, ranked by their log2FC**. It then asks whether genes belonging to a given pathway tend to cluster at the **top or bottom of the ranked list**. This makes GSEA:

- More sensitive (no arbitrary p-value threshold needed)
- Able to detect pathways where many genes shift slightly in the same direction

The output metric is the **Normalized Enrichment Score (NES)**:
- **NES > 0** — pathway genes are enriched at the top of the list (pathway is *activated* by dexamethasone)
- **NES < 0** — pathway genes are enriched at the bottom (pathway is *suppressed*)

#### Steps internally:

1. **Load all DEG results** — reads the full `DEG_results.csv` (not just significant)
2. **ID conversion** — same Ensembl → Entrez mapping as ORA
3. **Build ranked list** — named numeric vector `ENTREZID → log2FC`, sorted descending. Duplicate Entrez IDs kept by highest |LFC|.
4. **gsePathway()** — runs GSEA via the fgsea algorithm; `minGSSize=10`, `maxGSSize=500`, `eps=1e-10`
5. **Fallback** — if no pathways at padj < 0.05, retries at padj < 0.1
6. **Save results** + **generate plots**

#### GSEA output plots

| Plot | What it shows |
|------|--------------|
| `GSEA_ridgeplot` | Distribution of log2FC values for genes in each pathway; direction shows activation/suppression |
| `GSEA_dotplot` | Faceted (activated / suppressed) dotplot of top pathways |
| `GSEA_NES_barplot` | Horizontal bar chart of NES for top 10 activated + top 10 suppressed pathways |
| `GSEA_enrichment_1_*` ... `_4_*` | Running enrichment score plots for top 4 pathways — the classic GSEA "mountain plot" |
| `GSEA_emapplot` | Pathway similarity network coloured by NES |

---

### Step 9 — HTML report (RMarkdown)

**Module:** [`modules/report.nf`](modules/report.nf)  
**Script:** [`bin/report.R`](bin/report.R)  
**Template:** [`report/RNAseq_report.Rmd`](report/RNAseq_report.Rmd)  
**Tool:** rmarkdown, knitr, pandoc, DT  

Gathers all plots and result tables from DESeq2, Reactome ORA, and GSEA and renders them into a single **self-contained HTML file** (`results/report/RNAseq_analysis_report.html`).

The report includes:
- Interactive searchable DEG tables (DT package)
- All publication-quality plots embedded inline
- Session information and software versions
- Interpretation captions

You can share this single `.html` file with collaborators — it requires no internet connection to view.

---

### Step 10 — Quality summary (MultiQC)

**Module:** [`modules/multiqc.nf`](modules/multiqc.nf)  
**Tool:** MultiQC v1.21  

MultiQC collects the output logs from every step and aggregates them into a single interactive HTML dashboard. It automatically parses:

| Source | What is parsed |
|--------|---------------|
| FastQC | Per-sample quality metrics |
| fastp | Trimming statistics, before/after read counts |
| STAR logs | Mapping rates, read counts, splice junction stats |
| featureCounts `.summary` | Assigned vs. unassigned read counts |
| samtools flagstat | Mapped, paired, properly paired |
| samtools stats | Insert size distribution, base quality |

**Output:** `results/multiqc/multiqc_report.html`

> This is the first file you should open after a run. Any sample with unusual values immediately stands out in the colour-coded tables.

---

## Understanding every output file

```
results/
│
├── fastqc/                           ← Step 1 outputs
│   ├── SRR1039508_1_fastqc.html      ← Per-read-file QC report (R1)
│   ├── SRR1039508_1_fastqc.zip       ← Data for MultiQC
│   └── ...                           (2 files × 8 samples = 16 total)
│
├── fastp/                            ← Step 2 outputs
│   ├── SRR1039508.fastp.html         ← Trimming report per sample
│   ├── SRR1039508.fastp.json         ← Machine-readable trimming stats
│   └── ...                           (3 files × 8 samples = 24 total)
│
├── star/                             ← Step 3 outputs
│   └── SRR1039508/
│       ├── bam/
│       │   ├── SRR1039508.sortedByCoord.out.bam    ← Aligned reads
│       │   └── SRR1039508.sortedByCoord.out.bam.bai ← Index
│       ├── log/
│       │   └── SRR1039508.Log.final.out             ← Mapping statistics
│       └── spliceJunctions/
│           └── SRR1039508.SJ.out.tab                ← All detected splice junctions
│
├── bam_qc/                           ← Step 4 outputs
│   ├── SRR1039508.flagstat            ← Read counts: mapped, paired, etc.
│   ├── SRR1039508.stats               ← Detailed per-base quality stats
│   └── SRR1039508.idxstats            ← Per-chromosome read distribution
│
├── featurecounts/                    ← Step 5 outputs
│   ├── gene_counts.txt               ← THE COUNT MATRIX (input to DESeq2)
│   │                                    Columns: Geneid Chr Start End Strand Length
│   │                                             SRR1039508 SRR1039509 ... SRR1039521
│   └── gene_counts.txt.summary       ← Assignment stats per sample
│                                        Look for "Assigned" > 60%
│
├── deseq2/                           ← Step 6 outputs
│   ├── DEG_results.csv               ← ALL tested genes (sorted by padj)
│   │     Columns: gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
│   │     • baseMean: average normalised count across all samples
│   │     • log2FoldChange: treated vs control (apeglm-shrunk)
│   │     • lfcSE: standard error of log2FC
│   │     • stat: Wald test statistic
│   │     • padj: Benjamini-Hochberg adjusted p-value
│   │
│   ├── DEG_significant.csv           ← Filtered: padj<0.05 AND |log2FC|>1
│   │                                    Use this list for ORA and follow-up experiments
│   │
│   ├── normalized_counts.csv         ← Size-factor normalised counts
│   │                                    Use for comparing expression between samples
│   │
│   ├── vst_counts.csv                ← VST-transformed matrix
│   │                                    Use for PCA, clustering, heatmaps
│   │
│   ├── sessionInfo_deseq2.txt        ← R package versions (reproducibility)
│   │
│   ├── plots/
│   │   ├── PCA_plot.pdf/.png         ← Principal component analysis
│   │   ├── MA_plot.pdf/.png          ← Mean vs. fold change
│   │   ├── Volcano_plot.pdf/.png     ← Significance vs. fold change
│   │   ├── Heatmap_top_DEGs.pdf/.png ← Top 50 DEGs, Z-scored
│   │   ├── Sample_correlation_heatmap.pdf/.png ← Inter-sample distances
│   │   ├── Top9_DEG_counts.pdf/.png  ← Count distributions for top genes
│   │   └── Dispersion_plot.png       ← DESeq2 dispersion estimates
│   │
│   └── rds/
│       ├── dds.rds    ← DESeqDataSet — re-use in R for custom analyses
│       ├── res.rds    ← DESeqResults — the full results object
│       └── vsd.rds    ← VST SummarizedExperiment
│
├── pathway_analysis/
│   │
│   ├── reactome_ora/                 ← Step 7 outputs
│   │   ├── reactome_ORA_results.csv     ← All pathways for all DEGs
│   │   │     Columns: ID, Description, GeneRatio, BgRatio, pvalue,
│   │   │              p.adjust, qvalue, geneID, Count
│   │   │     • GeneRatio: DEGs in pathway / total DEGs
│   │   │     • BgRatio: all genes in pathway / all tested genes
│   │   │     • geneID: gene symbols of DEGs in this pathway
│   │   │
│   │   ├── reactome_ORA_upregulated.csv  ← Pathways from up-regulated DEGs only
│   │   ├── reactome_ORA_downregulated.csv ← Pathways from down-regulated DEGs only
│   │   ├── sessionInfo_reactome.txt
│   │   └── plots/
│   │       ├── ORA_dotplot.pdf/.png
│   │       ├── ORA_barplot.pdf/.png
│   │       ├── ORA_cnetplot.pdf/.png
│   │       ├── ORA_emapplot.pdf/.png
│   │       ├── ORA_upsetplot.pdf/.png
│   │       └── ORA_up_vs_down_dotplot.pdf/.png
│   │
│   └── reactome_gsea/               ← Step 8 outputs
│       ├── reactome_GSEA_results.csv    ← Enriched pathways (sorted by p.adjust)
│       │     Columns: ID, Description, setSize, enrichmentScore, NES,
│       │              pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment
│       │     • NES: Normalized Enrichment Score (>0 activated, <0 suppressed)
│       │     • rank: position in the ranked list where enrichment peaks
│       │     • leading_edge: tags=%, list=%, signal=% statistics
│       │     • core_enrichment: Entrez IDs of the "leading edge" genes
│       │
│       ├── sessionInfo_gsea.txt
│       └── plots/
│           ├── GSEA_ridgeplot.pdf/.png
│           ├── GSEA_dotplot.pdf/.png
│           ├── GSEA_NES_barplot.pdf/.png
│           ├── GSEA_enrichment_1_*.pdf/.png  ← Running score, pathway 1
│           ├── GSEA_enrichment_2_*.pdf/.png  ← Running score, pathway 2
│           ├── GSEA_enrichment_3_*.pdf/.png  ← Running score, pathway 3
│           ├── GSEA_enrichment_4_*.pdf/.png  ← Running score, pathway 4
│           └── GSEA_emapplot.pdf/.png
│
├── multiqc/                          ← Step 10 output
│   └── multiqc_report.html           ← THE QC DASHBOARD — open this first
│
├── report/                           ← Step 9 output
│   └── RNAseq_analysis_report.html   ← THE ANALYSIS REPORT — shareable
│
└── pipeline_info/                    ← Nextflow metadata
    ├── execution_timeline_*.html     ← Gantt chart of process execution times
    ├── execution_report_*.html       ← Resource usage (RAM, CPU, wall time)
    ├── execution_trace_*.txt         ← Tab-delimited log of every task
    └── pipeline_dag_*.svg            ← Directed acyclic graph of the workflow
```

---

## Repository structure

```
rna-dexamethasone/
├── main.nf                         # Main workflow entry point
├── nextflow.config                 # Pipeline configuration + profiles
│
├── modules/                        # DSL2 process modules (one tool = one file)
│   ├── fastqc.nf                   # FastQC quality assessment
│   ├── fastp.nf                    # fastp adapter trimming
│   ├── star_align.nf               # STAR splice-aware alignment
│   ├── bam_qc.nf                   # samtools flagstat / stats / idxstats
│   ├── featurecounts.nf            # featureCounts gene quantification
│   ├── deseq2.nf                   # DESeq2 module wrapper
│   ├── reactome.nf                 # Reactome ORA module wrapper
│   ├── gsea.nf                     # Reactome GSEA module wrapper
│   ├── multiqc.nf                  # MultiQC aggregation
│   └── report.nf                   # RMarkdown HTML report
│
├── subworkflows/                   # Logical pipeline segments
│   ├── qc_trim.nf                  # Calls FASTQC + FASTP
│   ├── align_quantify.nf           # Calls STAR + BAM_QC + FEATURECOUNTS
│   └── differential_expression.nf  # Calls DESEQ2 + REACTOME_ORA + REACTOME_GSEA + RENDER_REPORT
│
├── bin/                            # R analysis scripts called by modules
│   ├── deseq2.R                    # Full DESeq2 pipeline with all plots (414 lines)
│   ├── reactome.R                  # Reactome ORA (6 plot types)
│   ├── gsea.R                      # Reactome GSEA (5 plot types)
│   ├── report.R                    # RMarkdown rendering wrapper
│   ├── download_data.sh            # Download GSE52778 FASTQs (Linux)
│   └── download_data.ps1           # Download GSE52778 FASTQs (Windows/PowerShell)
│
├── report/
│   └── RNAseq_report.Rmd           # Comprehensive HTML report template
│
├── assets/
│   ├── samplesheet.csv             # Input sample manifest (8 GSE52778 samples)
│   ├── metadata.csv                # Sample metadata: condition, cell_line, batch, sex
│   └── multiqc_config.yml          # MultiQC branding and module config
│
├── conf/
│   └── envs/                       # Per-label conda environment definitions
│       ├── qc.yml                  # FastQC, fastp, MultiQC
│       ├── align.yml               # STAR, samtools, featureCounts
│       └── r_analysis.yml          # R + all Bioconductor packages
│
├── containers/
│   ├── Dockerfile                  # Builds rna-dexamethasone:1.0.0
│   └── environment.yml             # Full conda environment (alternative to Docker)
│
└── docs/
    ├── usage.md                    # Detailed usage guide
    └── pipeline_explained.md       # Deep-dive: every script, every line explained
```

---

## Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥ 23.04.0 | Workflow orchestration |
| Docker | ≥ 20.10 | Reproducible containerised execution |
| RAM | ≥ 32 GB | STAR genome loading (31 GB for GRCh38) |
| CPU | ≥ 8 cores | Parallel alignment |
| Disk | ≥ 150 GB | FASTQs (≈30 GB) + genome index (≈30 GB) + results |

> **WSL2 users (Windows):** The pipeline was developed and fully tested under WSL2 on Windows 11. Ensure your Docker Desktop is configured to use the WSL2 backend and that your data drive (e.g., `D:`) is mounted in WSL at `/mnt/d`.

---

## Installation

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version   # should print 23.x or higher
```

### 2. Clone the repository

```bash
git clone https://github.com/yourname/rna-dexamethasone.git
cd rna-dexamethasone
```

### 3. Build the Docker image

```bash
# This takes ~15–30 minutes (downloads all R packages)
docker build -t rna-dexamethasone:1.0.0 containers/

# Verify
docker run --rm rna-dexamethasone:1.0.0 STAR --version
docker run --rm rna-dexamethasone:1.0.0 Rscript -e "packageVersion('DESeq2')"
```

### 4. Download the GENCODE reference genome

```bash
mkdir -p refs

# GRCh38 primary assembly + GENCODE v46 annotation (~2 GB download)
wget -P refs/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
wget -P refs/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

gunzip refs/*.gz
```

### 5. Build the STAR genome index

```bash
# ⚠ Requires ~32 GB RAM and ~30 minutes
mkdir -p refs/star_index

STAR \
    --runMode genomeGenerate \
    --genomeDir refs/star_index \
    --genomeFastaFiles refs/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile refs/gencode.v46.annotation.gtf \
    --runThreadN 16 \
    --genomeSAindexNbases 14
```

### 6. Download the GSE52778 FASTQ files

**Linux/WSL (using SRA Toolkit):**
```bash
mkdir -p data

for SRR in SRR1039508 SRR1039509 SRR1039512 SRR1039513 \
            SRR1039516 SRR1039517 SRR1039520 SRR1039521; do
    echo "Downloading $SRR..."
    fasterq-dump --split-files --gzip --threads 8 --outdir data/ $SRR
done
```

**Windows PowerShell (using the included script):**
```powershell
.\bin\download_data.ps1
```

---

## Running the pipeline

### With Docker (recommended)

```bash
nextflow run main.nf \
    -profile docker \
    --star_index refs/star_index \
    --gtf refs/gencode.v46.annotation.gtf \
    --outdir results \
    -resume
```

> `-resume` tells Nextflow to reuse cached results from previous runs. Add this flag every time — it makes re-runs after any failure very fast.

### With Conda

```bash
nextflow run main.nf \
    -profile conda \
    --star_index refs/star_index \
    --gtf refs/gencode.v46.annotation.gtf \
    --outdir results
```

### On SLURM HPC

```bash
nextflow run main.nf \
    -profile slurm,docker \
    --star_index /shared/refs/GRCh38/star_index \
    --gtf /shared/refs/GRCh38/gencode.v46.annotation.gtf \
    --outdir /scratch/$USER/rna-dex-results \
    -resume \
    -bg   # run in background
```

### Run R scripts standalone (skip alignment)

If you already have a count matrix (e.g., from the `airway` R package or another pipeline):

```bash
# Step 1: DESeq2
Rscript bin/deseq2.R \
    --counts  results/featurecounts/gene_counts.txt \
    --metadata assets/metadata.csv \
    --pvalue 0.05 --lfc 1.0 \
    --outdir results/deseq2

# Step 2: Reactome ORA (requires DEG_significant.csv from step 1)
Rscript bin/reactome.R \
    --deg results/deseq2/DEG_significant.csv \
    --organism human \
    --pvalue 0.05 \
    --outdir results/pathway_analysis/reactome_ora

# Step 3: Reactome GSEA (requires DEG_results.csv — the full table)
Rscript bin/gsea.R \
    --deg results/deseq2/DEG_results.csv \
    --organism human \
    --pvalue 0.05 \
    --outdir results/pathway_analysis/reactome_gsea
```

---

## Parameters reference

### Input / Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `assets/samplesheet.csv` | CSV with columns: `sample, condition, fastq1, fastq2` |
| `--metadata` | `assets/metadata.csv` | Sample metadata for DESeq2 `colData` |
| `--outdir` | `results` | Output directory |

### Reference genome

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--star_index` | `null` | Path to pre-built STAR genome directory |
| `--gtf` | `null` | Genome annotation GTF file (GENCODE recommended) |
| `--genome` | `GRCh38` | Genome assembly name (used for logging only) |

### Analysis thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pvalue_cutoff` | `0.05` | BH-adjusted p-value threshold for DEGs and enrichment |
| `--lfc_threshold` | `1.0` | Minimum \|log2FC\| for significant DEGs (= 2-fold change) |
| `--organism` | `human` | Organism for ReactomePA (`human` or `mouse`) |
| `--deseq2_min_count` | `10` | Minimum total count per gene across all samples |
| `--deseq2_lfc_shrink` | `true` | Apply apeglm LFC shrinkage (strongly recommended) |

### Tool-specific

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastp_min_length` | `30` | Minimum read length after trimming |
| `--fastp_quality` | `20` | Phred quality threshold for trimming |
| `--star_threads` | `4` | Threads for STAR alignment |
| `--star_two_pass` | `false` | Enable STAR 2-pass mode (needs more RAM) |
| `--fc_threads` | `4` | Threads for featureCounts |
| `--fc_feature_type` | `exon` | GTF feature to count over |
| `--fc_group_attribute` | `gene_id` | GTF attribute to summarise by |

### Resource limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Hard cap on CPUs per process |
| `--max_memory` | `12.GB` | Hard cap on memory per process |
| `--max_time` | `240.h` | Hard cap on wall-clock time per process |

---

## Samplesheet format

```csv
sample,condition,fastq1,fastq2
SRR1039508,control,data/SRR1039508_1.fastq.gz,data/SRR1039508_2.fastq.gz
SRR1039509,treated,data/SRR1039509_1.fastq.gz,data/SRR1039509_2.fastq.gz
```

**Rules:**
- `sample` — unique identifier; becomes the BAM file name and column name in the count matrix
- `condition` — must include a level called exactly `control` (this becomes the DESeq2 reference level)
- `fastq1` — path to R1; can be absolute or relative to the project directory
- `fastq2` — path to R2; leave empty for single-end data

---

## Expected results (GSE52778)

Based on the original publication (Himes et al., 2014) and the full pipeline run on this dataset:

| Metric | Expected value |
|--------|----------------|
| Mapping rate | > 90% per sample |
| Genes tested | ~22,000 after pre-filtering |
| Significant DEGs (padj < 0.05, \|LFC\| > 1) | ~4,000–6,000 |
| **Top up-regulated genes** | FKBP5, CRISPLD2, DUSP1, ZBTB16, TSC22D3 |
| **Top down-regulated genes** | IL8 (CXCL8), CXCL1, CCL2, TNF, IL6 |
| **Key activated pathways** | Glucocorticoid receptor signalling, Gene expression (transcription) |
| **Key suppressed pathways** | Interleukin signalling, Innate immune system, Cytokine signalling |

### Biological interpretation

The results tell a consistent story:

1. **FKBP5, DUSP1, TSC22D3** — canonical GR target genes — are strongly up-regulated. This confirms successful GR activation.
2. **IL8, CXCL1, TNF, CCL2** — pro-inflammatory cytokines — are suppressed. This is the mechanism of dexamethasone's anti-inflammatory effect.
3. **Reactome ORA** confirms suppression of the "Interleukin-1 signalling" and "MAPK cascade" pathways.
4. **GSEA** shows the full transcriptome-level picture: the entire glucocorticoid signalling gene set is shifted upward (NES > 0), and cytokine signalling gene sets are shifted downward (NES < 0).

---

## Troubleshooting

**`No matching samples between count matrix and metadata`**
> The sample names in the count matrix header must match the `sample` column in `metadata.csv`. The script strips `.sortedByCoord.out.bam` and path prefixes automatically. Check that your SRR IDs are consistent.

**STAR runs out of memory**
> STAR requires ~31 GB to load the GRCh38 genome. Increase `--max_memory '64.GB'` or, if on a laptop, reduce genome chromosomes. On WSL2, increase the memory limit in `.wslconfig` under `[wsl2]` → `memory=32GB`.

**`enrichPathway` finds 0 pathways**
> Relax thresholds: `--pvalue_cutoff 0.1` or `--lfc_threshold 0.5`. If Ensembl ID → Entrez mapping fails, check that gene IDs in your count matrix are in `ENSG...` format (GENCODE) rather than gene symbols.

**GSEA: `Insufficient number of permutations`**
> GSEA needs a minimum number of genes to be reliable. Use the full `DEG_results.csv` (all tested genes), not the filtered significant table.

**featureCounts crashes with `Segmentation fault`**
> This can occur when the Docker-compiled binary and the host kernel version are incompatible (WSL2-specific). The module is configured to use a host-native binary (`/home/yer_kanat/Downloads/subread-2.0.6-source/bin/featureCounts`). See [docs/usage.md](docs/usage.md) for instructions on compiling featureCounts from source.

**Docker: `Error response from daemon: cannot start a stopped process`**
> Restart Docker Desktop and ensure the WSL2 backend is enabled.

---

## Citations

If you use this pipeline in your research, please cite:

- **DESeq2:** Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology. 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8
- **ReactomePA:** Yu G, He QY. *ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization.* Mol BioSyst. 2016;12(2):477-479. https://doi.org/10.1039/C5MB00663E
- **clusterProfiler:** Wu T, et al. *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.* Innovation. 2021;2(3):100141. https://doi.org/10.1016/j.xinn.2021.100141
- **STAR:** Dobin A, et al. *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics. 2013;29(1):15-21. https://doi.org/10.1093/bioinformatics/bts635
- **featureCounts:** Liao Y, et al. *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.* Bioinformatics. 2014;30(7):923-930. https://doi.org/10.1093/bioinformatics/btt656
- **fastp:** Chen S, et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics. 2018;34(17):i884-i890. https://doi.org/10.1093/bioinformatics/bty560
- **MultiQC:** Ewels P, et al. *MultiQC: summarize analysis results for multiple tools and samples in a single report.* Bioinformatics. 2016;32(19):3047-3048. https://doi.org/10.1093/bioinformatics/btw354
- **GSE52778 dataset:** Himes BE, et al. *RNA-seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells.* PLOS ONE. 2014;9(6):e99625. https://doi.org/10.1371/journal.pone.0099625
- **Nextflow:** Di Tommaso P, et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology. 2017;35:316-319. https://doi.org/10.1038/nbt.3820

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

<div align="center">
Built with ❤️ using Nextflow DSL2, R/Bioconductor, and open-source bioinformatics tools.<br>
Data: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778">GSE52778</a> (Himes et al., PLOS ONE, 2014)
</div>
