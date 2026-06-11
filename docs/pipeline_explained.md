# Pipeline Deep-Dive: Everything Explained

This document explains **every file, every script, every parameter, and every output** of the `rna-dexamethasone` Nextflow pipeline. It is written for someone who understands basic molecular biology but may be new to bioinformatics software.

---

## Table of Contents

1. [How Nextflow works](#how-nextflow-works)
2. [main.nf — the entry point](#mainnf--the-entry-point)
3. [nextflow.config — the control panel](#nextflowconfig--the-control-panel)
4. [The Docker container](#the-docker-container)
5. [Subworkflow: QC_TRIM](#subworkflow-qc_trim)
   - [modules/fastqc.nf](#modulesfastqcnf)
   - [modules/fastp.nf](#modulesfaspnf)
6. [Subworkflow: ALIGN_QUANTIFY](#subworkflow-align_quantify)
   - [modules/star_align.nf](#modulesstar_alignnf)
   - [modules/bam_qc.nf](#modulesbam_qcnf)
   - [modules/featurecounts.nf](#modulesfeaturecountsnf)
7. [Subworkflow: DIFF_EXPRESSION](#subworkflow-diff_expression)
   - [bin/deseq2.R](#bindeseq2r)
   - [bin/reactome.R](#binreactomer)
   - [bin/gsea.R](#bingseaR)
   - [bin/report.R](#binreportr)
8. [modules/multiqc.nf](#modulesmultiqcnf)
9. [Assets: samplesheet.csv and metadata.csv](#assets-samplesheetcsv-and-metadatacsv)
10. [How to interpret each output file](#how-to-interpret-each-output-file)

---

## How Nextflow works

Nextflow is a **workflow manager**. You describe a pipeline as a set of **processes** (each runs a shell command or script), connected by **channels** (queues of data). Nextflow then figures out which processes can run in parallel, submits them to the executor (your local machine or an HPC cluster), and handles caching, retries, and logging automatically.

Key concepts:

| Concept | What it is |
|---------|-----------|
| **Process** | One computational step. Has `input`, `output`, and `script` blocks. |
| **Channel** | A queue of items (files, strings, tuples) that flow between processes. |
| **Subworkflow** | A named group of processes that together form a logical stage. |
| **Profile** | A named set of configuration overrides (e.g., `docker`, `conda`, `slurm`). |
| **-resume** | Nextflow fingerprints each task's inputs. If inputs have not changed, the cached result is reused. |
| **DSL2** | The current Nextflow syntax (version 2). Processes can be imported as modules and reused. |

When you run `nextflow run main.nf -profile docker`, Nextflow:
1. Reads `nextflow.config` and applies the `docker` profile
2. Parses `main.nf` to build the workflow DAG
3. For each process, starts a Docker container with `rna-dexamethasone:1.0.0`
4. Copies input files into the container's work directory
5. Runs the script inside the container
6. Copies declared outputs back to `results/`

---

## main.nf — the entry point

```
main.nf
```

This 123-line file is the top of the pipeline. Here is what each section does:

### Imports

```groovy
include { QC_TRIM          } from './subworkflows/qc_trim'
include { ALIGN_QUANTIFY   } from './subworkflows/align_quantify'
include { DIFF_EXPRESSION  } from './subworkflows/differential_expression'
include { MULTIQC          } from './modules/multiqc'
```

DSL2 `include` statements import named workflows and processes from other files. This keeps `main.nf` clean and each logical step in its own file.

### Parameter validation

```groovy
def validateParams() {
    if (!params.samplesheet) error "..."
    if (!params.star_index && !params.genome) error "..."
    if (!params.gtf) error "..."
}
```

Fails immediately with a human-readable error if required parameters are missing, before any compute is wasted.

### Input channel

```groovy
Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true, strip: true)
    .map { row ->
        def meta = [ id: row.sample, condition: row.condition,
                     single_end: row.fastq2 ? false : true ]
        def r1 = file(row.fastq1, checkIfExists: true)
        def r2 = row.fastq2 ? file(row.fastq2, checkIfExists: true) : []
        return [meta, r1, r2]
    }
    .set { ch_reads }
```

This reads `samplesheet.csv` row-by-row and converts each row into a **tuple** `[meta, r1, r2]` where `meta` is a Groovy map of sample attributes. The `ch_reads` channel emits one tuple per sample. All 8 samples are emitted in parallel into downstream processes.

### MultiQC channel construction

```groovy
ch_multiqc_files = Channel.empty()
    .mix(QC_TRIM.out.fastqc_zip.map{ meta, files -> files }.collect())
    .mix(QC_TRIM.out.trim_log.collect())
    .mix(ALIGN_QUANTIFY.out.star_log.map{ meta, file -> file }.collect())
    ...
    .collect()
```

`mix()` merges multiple channels. `collect()` waits until all items arrive and bundles them into a single list. This is how all QC logs from all steps and all samples get gathered into one place for MultiQC to process.

---

## nextflow.config — the control panel

This file configures every aspect of pipeline execution. Key sections:

### params block

Defines all default parameter values. Any of these can be overridden on the command line with `--param_name value`.

```groovy
params {
    samplesheet    = "${projectDir}/assets/samplesheet.csv"  // default input
    pvalue_cutoff  = 0.05                                     // DESeq2 + enrichment
    lfc_threshold  = 1.0                                      // 2-fold minimum
    max_memory     = '12.GB'                                  // per-process cap
    ...
}
```

### check_max function

```groovy
check_max = { obj, type ->
    if (type == 'memory') {
        if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
            return params.max_memory as nextflow.util.MemoryUnit
        else return obj
    }
```

This closure is used throughout as `params.check_max(64.GB, 'memory')`. It returns whichever is smaller: the requested amount or the `max_memory` cap. This prevents processes from requesting more than your machine has.

### process block — resource labels

```groovy
withLabel: 'process_high' { cpus = 4; memory = '10.GB'; time = '12.h' }
withLabel: 'process_r'    { cpus = 2; memory = '8.GB';  time = '4.h'  }
```

Each module declares a `label` (e.g., `label 'process_high'`). The config matches labels to resource allocations. R analysis uses `process_r` because it needs plenty of RAM for loading Bioconductor data packages.

### Per-process overrides

```groovy
withName: 'STAR_ALIGN' {
    memory = { params.check_max(64.GB * task.attempt, 'memory') }
}
```

The `task.attempt` multiplier means: if STAR fails (OOM), retry with double the memory.

### Docker block

```groovy
docker {
    enabled       = true
    userEmulation = true
    runOptions    = '-u $(id -u):$(id -g) -v /mnt/d:/mnt/d'
}
```

`userEmulation = true` runs the container as your host user, preventing root-owned output files. The `-v /mnt/d:/mnt/d` bind mount exposes the D: drive inside the container (WSL2 setup).

### Profiles

| Profile | When to use |
|---------|------------|
| `docker` | Local machine or cloud with Docker installed |
| `conda` | When Docker is unavailable; uses mamba to install environments |
| `singularity` | HPC environments that use Singularity instead of Docker |
| `slurm` | HPC clusters with the SLURM job scheduler |
| `test` | Quick test run with a small dataset |

### Reporting config

```groovy
timeline { file = "${params.outdir}/pipeline_info/execution_timeline_*.html" }
trace    { file = "...execution_trace_*.txt" }
dag      { file = "...pipeline_dag_*.svg" }
```

Nextflow automatically generates:
- **Timeline** — a Gantt chart showing which processes ran when and how long they took
- **Trace** — a TSV with every task's CPU/RAM usage
- **DAG** — the workflow directed acyclic graph (dependency diagram)
- **Report** — resource usage charts

---

## The Docker container

**File:** `containers/Dockerfile`

Built from `rocker/r-ver:4.5` (a Debian Linux image with R 4.5 pre-installed). The build proceeds in layers:

### Layer 1 — System libraries
```dockerfile
RUN apt-get install -y wget curl samtools openjdk-17-jre-headless \
    libcurl4-openssl-dev libxml2-dev ...
```
- `samtools` — BAM file manipulation
- `openjdk-17` — required by FastQC (Java application)
- `libcurl4-openssl-dev`, `libssl-dev`, `libxml2-dev` — required to compile R packages that make HTTP requests (BiocManager, biomaRt)
- `libharfbuzz-dev`, `libpng-dev`, etc. — required for R graphics output (ggplot2 → ragg → PNG)

### Layer 2 — FastQC
```dockerfile
ARG FASTQC_VERSION=0.12.1
RUN wget https://...fastqc_v0.12.1.zip && unzip && chmod +x && ln -s ...
```
Downloads the FastQC Java wrapper, makes it executable, and symlinks it to `/usr/local/bin/fastqc`.

### Layer 3 — fastp
```dockerfile
ARG FASTP_VERSION=0.23.4
RUN wget https://...fastp -O /usr/local/bin/fastp && chmod +x
```
fastp is distributed as a single static binary — the simplest install possible.

### Layer 4 — STAR
```dockerfile
ARG STAR_VERSION=2.7.11b
RUN wget https://...STAR_2.7.11b.zip && unzip && cp Linux_x86_64_static/STAR /usr/local/bin/
```
Uses the pre-compiled static binary from the STAR GitHub releases. No compilation needed.

### Layer 5 — featureCounts (compiled from source)
```dockerfile
ARG SUBREAD_VERSION=2.0.6
RUN wget .../subread-2.0.6-source.tar.gz \
    && tar xzf ... \
    && cd subread-2.0.6-source/src \
    && make -f Makefile.Linux
```
Compiled from source. This is the most reliable approach on WSL2 where pre-built binaries can segfault due to kernel version mismatches.

### Layer 6 — MultiQC
```dockerfile
RUN pip3 install --no-cache-dir multiqc==1.21
```
MultiQC is a Python package, installed via pip.

### Layer 7 — R CRAN packages
```dockerfile
RUN R -e "install.packages(c('optparse','ggplot2','dplyr','pheatmap',
    'rmarkdown','knitr','DT','cowplot','ggrepel','stringr','viridis'), ...)"
```
Key packages:
- `optparse` — command-line argument parsing for the R scripts
- `ggplot2` — all plotting
- `pheatmap` — heatmaps
- `rmarkdown` + `knitr` — HTML report generation
- `DT` — interactive DataTables in HTML
- `ggrepel` — non-overlapping labels in volcano/PCA plots

### Layer 8 — Bioconductor packages
```dockerfile
RUN R -e "BiocManager::install(c('DESeq2','ReactomePA','clusterProfiler',
    'enrichplot','org.Hs.eg.db','EnhancedVolcano','fgsea','ggupset'), ...)"
```
Key packages:
- `DESeq2` — the differential expression engine
- `ReactomePA` — Reactome pathway enrichment (ORA and GSEA)
- `clusterProfiler` — the backbone of enrichment analysis
- `enrichplot` — all pathway visualisations (dotplot, cnetplot, emapplot, ridgeplot)
- `org.Hs.eg.db` — human gene ID annotation database (Ensembl → Entrez → Symbol)
- `EnhancedVolcano` — publication-quality volcano plots
- `fgsea` — fast GSEA implementation used by gsePathway
- `ggupset` — upset plots

### Layer 9 — pandoc
```dockerfile
RUN apt-get install -y perl pandoc
```
Added at the end (separate layer) to avoid invalidating the long R package installation layers. `pandoc` is the document converter that `rmarkdown::render()` uses to convert the knitted HTML from Markdown.

---

## Subworkflow: QC_TRIM

**File:** `subworkflows/qc_trim.nf`

Takes `ch_reads` (8 tuples of `[meta, r1, r2]`) and runs FastQC + fastp on each in parallel.

```groovy
FASTQC(reads)      // runs on raw reads — captures pre-trim quality
FASTP(reads)       // trims adapters and low-quality ends
```

Both modules receive the **original raw reads** (not the trimmed ones). FastQC thus shows the quality of the data **before** trimming, which is important for QC diagnostics.

### modules/fastqc.nf

```groovy
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    input:  tuple val(meta), path(r1), path(r2)
    output: tuple val(meta), path("*.html"), emit: html
            tuple val(meta), path("*.zip"),  emit: zip
    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}
```

- `tag "${meta.id}"` — labels each task in the Nextflow log with the sample name (e.g., `SRR1039508`)
- `label 'process_low'` — request 2 CPUs / 4 GB RAM
- `path("*.html")` — Nextflow captures any `.html` file produced in the work directory
- FastQC automatically names outputs after the input file: `SRR1039508_1_fastqc.html`

### modules/fastp.nf

fastp receives paired-end reads and produces:
- `*.R1.trimmed.fastq.gz` and `*.R2.trimmed.fastq.gz` — the cleaned reads
- `*.fastp.json` / `*.fastp.html` — trimming statistics
- `*.fastp.log` — the log file parsed by MultiQC

Key fastp flags (from `nextflow.config`):
- `--qualified_quality_phred 20` — trim bases with Q < 20
- `--length_required 30` — discard reads shorter than 30 bp after trimming
- `-w ${task.cpus}` — multi-threaded

---

## Subworkflow: ALIGN_QUANTIFY

**File:** `subworkflows/align_quantify.nf`

Receives 8 trimmed FASTQ pairs from `QC_TRIM.out.trimmed_reads` and runs alignment + QC + counting.

### modules/star_align.nf

STAR runs once per sample, inside Docker. The critical STAR parameters:

| Flag | Value | Purpose |
|------|-------|---------|
| `--genomeDir` | `params.star_index` | Pre-built genome index directory |
| `--readFilesCommand zcat` | — | Input FASTQs are gzipped; `zcat` decompresses on the fly |
| `--outSAMtype BAM SortedByCoordinate` | — | Write coordinate-sorted BAM directly (no separate sort step) |
| `--outSAMattributes NH HI AS NM MD` | — | NH=multi-mapping count, HI=hit index, AS=alignment score, NM=mismatches, MD=mismatch positions |
| `--outSAMstrandField intronMotif` | — | Adds strand info for stranded libraries |
| `--outFilterIntronMotifs RemoveNoncanonical` | — | Drop reads with unusual splice junctions |
| `--limitBAMsortRAM 1500000000` | 1.5 GB | Cap sorting memory to prevent WSL OOM |

After STAR, `samtools index` creates the `.bai` file. The BAM and BAI are renamed to strip STAR's verbose `Aligned.` prefix.

**Log.final.out — how to read it:**
```
                          Number of input reads | 25000000
                      Average input read length | 150
                    Uniquely mapped reads number | 23500000
                         Uniquely mapped reads % | 94.00%
              Number of reads mapped to multiple loci | 800000
                   % of reads mapped to multiple loci | 3.20%
```
- Aim for > 90% uniquely mapped
- > 5% multi-mappers can indicate repetitive regions or rRNA contamination
- > 5% unmapped usually means wrong reference or heavily degraded RNA

### modules/bam_qc.nf

Runs three samtools commands in one process per sample:

```bash
samtools flagstat ${bam} > ${prefix}.flagstat
samtools stats    ${bam} > ${prefix}.stats
samtools idxstats ${bam} > ${prefix}.idxstats
```

**Flagstat output example:**
```
24300000 + 0 in total (QC-passed reads + QC-failed reads)
24300000 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
23500000 + 0 mapped (96.71% : N/A)
24300000 + 0 paired in sequencing
12150000 + 0 read1
12150000 + 0 read2
23000000 + 0 properly paired (94.65% : N/A)
```
- "properly paired" means both R1 and R2 mapped in the expected orientation and distance — the most informative metric

**Idxstats — example:**
```
chr1    248956422   4500000   0
chr2    242193529   3900000   0
...
chrM    16569        12000   0    ← mitochondrial
```
If `chrM` has a disproportionate fraction (> 5%), the sample may have high mitochondrial RNA contamination.

### modules/featurecounts.nf

featureCounts is the **only process not running inside Docker** (it uses a host-compiled binary at a hardcoded path). This was a deliberate workaround for WSL2 segfaults with the container-compiled binary.

```groovy
input:
val bam_files   // NOTE: 'val' not 'path' — receives pre-resolved absolute path strings
val gtf         // absolute path to GTF file
```

The `val` input type tells Nextflow not to stage-copy the BAM files into the work directory — featureCounts accesses them directly via their absolute paths. This matters because the BAMs are large (>5 GB each) and copying all 8 would take hours.

**featureCounts key flags:**
- `-t exon` — count reads that overlap **exon** features in the GTF (not gene body or intron)
- `-g gene_id` — after counting per-exon, summarise to gene level using the `gene_id` attribute
- `-B` — paired-end: both reads in a pair must map; a read pair counts once
- `-C` — discard chimeric pairs (R1 and R2 on different chromosomes)
- `--minOverlap 10` — a read must overlap a feature by ≥ 10 bp

**gene_counts.txt format:**
```
# Program:featureCounts v2.0.6; ...
Geneid         Chr  Start   End      Strand  Length  SRR1039508  SRR1039509  ...
ENSG00000000003 chrX 99891803 99894988 - 4535    850         923         ...
```
- First 6 columns: gene annotation (ignored by DESeq2)
- Remaining columns: raw count per sample

**gene_counts.txt.summary format:**
```
Status                     SRR1039508  SRR1039509 ...
Assigned                   20000000    20500000   ...
Unassigned_Unmapped        100000      95000      ...
Unassigned_MultiMapping    500000      480000     ...
Unassigned_NoFeatures      800000      790000     ...
Unassigned_Ambiguity       200000      210000     ...
```
`Assigned` should be > 60% for a clean sample. High `Unassigned_NoFeatures` suggests a mismatch between the GTF and genome versions, or intronic reads (expected in some cell types).

---

## Subworkflow: DIFF_EXPRESSION

**File:** `subworkflows/differential_expression.nf`

Receives the count matrix and metadata file, runs 4 sequential modules.

### bin/deseq2.R

This is the statistical heart of the pipeline. Here is what each section does:

#### Loading and cleaning (lines 45–78)

```r
counts_raw <- read.delim(opt$counts, comment.char="#", check.names=FALSE)
gene_info  <- counts_raw[, 1:6]     # annotation columns
count_mat  <- as.matrix(counts_raw[, 7:ncol(counts_raw)])  # count columns
```

featureCounts puts a `#` comment line at the top. `comment.char="#"` skips it. `check.names=FALSE` preserves column names with special characters (e.g., slashes in paths).

```r
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))      # remove path
colnames(count_mat) <- gsub("\\.bam$", "", colnames(count_mat))   # remove .bam
colnames(count_mat) <- gsub("\\.sortedByCoord\\.out$", "", ...)   # remove suffix
```

The count matrix column names look like `/mnt/d/star/SRR1039508.sortedByCoord.out.bam`. These three `gsub` calls strip it down to just `SRR1039508` to match the `sample` column in `metadata.csv`.

#### DESeqDataSet creation (lines 82–92)

```r
metadata$condition <- relevel(metadata$condition, ref="control")
dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = metadata,
    design    = ~ condition
)
```

`relevel(..., ref="control")` ensures DESeq2 computes **treated ÷ control** (positive LFC = up in treated). Without this, the direction could be reversed.

#### Pre-filtering (lines 96–102)

```r
keep <- rowSums(counts(dds)) >= opt$min_count   # default: >= 10
dds  <- dds[keep, ]
```

Removes genes with negligible total signal. This is not arbitrary — DESeq2's statistical model assumes reasonable counts, and including near-zero count genes inflates the multiple testing burden (Bonferroni/BH corrections are based on the total number of tests).

#### DESeq2 fitting (line 107)

```r
dds <- DESeq(dds, parallel=FALSE)
```

Internally, `DESeq()` performs three steps:
1. `estimateSizeFactors` — normalise for library size differences
2. `estimateDispersions` — fit gene-wise dispersion with shrinkage
3. `nbinomWaldTest` — Wald test: is log2FC ≠ 0?

#### LFC shrinkage (lines 124–135)

```r
coef_name <- resultsNames(dds)[2]   # e.g. "condition_treated_vs_control"
res <- lfcShrink(dds, coef=coef_name, type="apeglm")
```

apeglm is a Bayesian maximum a posteriori (MAP) estimator for log fold change. It places a Cauchy prior on the LFC, which:
- Shrinks LFCs for low-count genes towards zero (they have too little evidence)
- Leaves large LFCs for high-count genes nearly unchanged (they have strong evidence)

The result: the volcano plot shows a clean, interpretable picture rather than an explosion of extreme-LFC low-count genes.

#### Output tables

| File | Contents | Use for |
|------|----------|---------|
| `DEG_results.csv` | All ~22,000 genes, sorted by padj | GSEA (needs the full ranked list) |
| `DEG_significant.csv` | padj<0.05 AND \|LFC\|>1 | ORA, follow-up experiments |
| `normalized_counts.csv` | Size-factor normalised counts | Comparing expression between samples |
| `vst_counts.csv` | VST-transformed | PCA, heatmaps, clustering |
| `dds.rds`, `res.rds`, `vsd.rds` | R objects | Custom downstream analysis in R |

#### PCA plot — how to read it

Principal Component Analysis (PCA) is a dimensionality reduction method. It compresses thousands of gene dimensions into two orthogonal axes (PC1, PC2) that capture the most variance.

- **Ideal result:** The 8 samples form two tight clusters (control vs. treated) that are well separated on PC1. PC1 should explain > 50% of variance if the treatment effect dominates.
- **Watch out for:** Samples from the same donor clustering together on PC2 (batch effect by cell line). This is expected and normal — the paired design in DESeq2 handles it.

#### Volcano plot — how to read it

- X-axis: log2 fold change (shrunk). 0 = no change. +1 = 2-fold up. -1 = 2-fold down.
- Y-axis: -log10(adjusted p-value). Higher = more significant. 1.3 = padj 0.05.
- Red dots (right): significantly up-regulated (padj<0.05, LFC>1)
- Blue dots (left): significantly down-regulated (padj<0.05, LFC<-1)
- Grey: not significant
- Top 30 most significant genes are labelled

#### Heatmap — how to read it

- Rows: top 50 significant DEGs (by padj)
- Columns: 8 samples
- Values: Z-scores of VST counts (how many SDs above or below each gene's mean across all samples)
- Red = higher than mean; Blue = lower than mean
- The two condition clusters (control left, treated right) should be clearly separated

---

### bin/reactome.R

#### ID conversion — why it is needed

DESeq2 outputs Ensembl gene IDs (ENSG...) because that is what GENCODE GTF uses. ReactomePA internally uses **Entrez IDs** (integer IDs from NCBI). The `bitr()` function from clusterProfiler converts between them using the `org.Hs.eg.db` annotation package (a pre-built SQLite database).

```r
entrez_df <- bitr(gene_ids, fromType="ENSEMBL", toType="ENTREZID",
                  OrgDb=org.Hs.eg.db, drop=TRUE)
```

Some Ensembl IDs do not have an Entrez equivalent (novel genes, pseudogenes) — `drop=TRUE` silently removes them. The mapping rate is printed to the log.

#### Over-Representation Analysis — the maths

For each Reactome pathway P with k genes, given a DEG list of n genes from a universe of N total tested genes, the ORA p-value is the probability under the hypergeometric distribution:

```
P(X >= q) where X ~ Hypergeometric(N, K, n)
```
- N = total genes tested (~22,000)
- K = genes in pathway P
- n = number of DEGs
- q = DEGs that are also in pathway P

This is equivalent to Fisher's exact test for a 2×2 contingency table. BH correction is applied across all pathways tested.

#### Running ORA three times

```r
enrich_all  <- run_enrichPathway(entrez_ids, "All DEGs")
enrich_up   <- run_enrichPathway(up_entrez,  "Up-regulated")
enrich_down <- run_enrichPathway(down_entrez, "Down-regulated")
```

Three separate analyses:
1. **All DEGs** — general picture; used for most plots
2. **Up-regulated only** — which pathways are **activated** by dexamethasone?
3. **Down-regulated only** — which pathways are **suppressed**?

#### ORA result columns

| Column | Meaning |
|--------|---------|
| `ID` | Reactome pathway ID (e.g., R-HSA-5663213) |
| `Description` | Human-readable pathway name |
| `GeneRatio` | `q/n` — fraction of your DEGs in this pathway |
| `BgRatio` | `K/N` — fraction of all genes in this pathway |
| `pvalue` | Raw hypergeometric p-value |
| `p.adjust` | BH-adjusted p-value |
| `qvalue` | Storey q-value (alternative FDR estimate) |
| `geneID` | Gene symbols of your DEGs in this pathway (e.g., `FKBP5/DUSP1/IL8`) |
| `Count` | Number of DEGs in this pathway (= q) |

---

### bin/gsea.R

#### The ranked gene list

```r
gene_list <- setNames(deg_mapped$log2FoldChange, deg_mapped$ENTREZID)
gene_list <- sort(gene_list, decreasing=TRUE)
```

Every tested gene gets a score equal to its log2FC (from the full `DEG_results.csv`, not just significant genes). Sorting descending puts the most up-regulated genes first. This creates a ranking from "most activated by dexamethasone" to "most suppressed."

#### GSEA algorithm (fgsea)

For each Reactome pathway:
1. Mark which positions in the ranked list contain pathway genes
2. Walk down the list, adding to a running score when you hit a pathway gene and subtracting when you hit a non-pathway gene
3. The **Enrichment Score (ES)** is the maximum deviation of this running sum from zero
4. The **NES** normalises ES by the mean ES from permuted gene lists, making it comparable across pathways of different sizes

#### GSEA result columns

| Column | Meaning |
|--------|---------|
| `NES` | Normalized Enrichment Score. >0 = activated, <0 = suppressed |
| `rank` | Position in ranked list where ES peaks (the "leading edge point") |
| `leading_edge` | `tags=X%` (fraction of pathway genes in leading edge), `list=Y%` (list fraction), `signal=Z%` (combined) |
| `core_enrichment` | Entrez IDs of genes in the leading edge (the "driving" genes) |

#### Enrichment score plot — reading it

The classic GSEA mountain plot (saved as `GSEA_enrichment_1_*` etc.):
- Top panel: the running enrichment score along the ranked list
- Middle panel: vertical tick marks showing where each pathway gene falls in the ranked list
- Bottom panel: the ranked metric (log2FC) value at each position
- A mountain curving up and peaking at the left = pathway genes are at the top of the list = activated
- A valley curving down and reaching minimum at the right = pathway genes are at the bottom = suppressed

---

### bin/report.R

```r
rmarkdown::render(
    input       = rmd_path,
    output_file = output_html,
    output_dir  = opt$outdir,
    params      = list(deseq2_dir = ..., ora_dir = ..., gsea_dir = ...)
)
```

`rmarkdown::render()` processes the `.Rmd` template:
1. `knitr` executes each R code chunk embedded in the template
2. Each chunk reads and plots the pre-computed results
3. The resulting Markdown (with embedded base64-encoded images) is passed to `pandoc`
4. `pandoc` converts it to a self-contained HTML file

The `output-format: self_contained: true` option encodes all images and CSS inline — the resulting `.html` file is completely standalone.

---

## modules/multiqc.nf

MultiQC runs **once**, after all other processes complete (enforced by the `collect()` call in `main.nf` which waits for all files before emitting).

```bash
multiqc . \
    --filename multiqc_report.html \
    --outdir ${params.outdir}/multiqc \
    --config ${projectDir}/assets/multiqc_config.yml
```

MultiQC auto-discovers all supported log files in the current directory and its subdirectories. It knows how to parse > 100 bioinformatics tools, including every tool used in this pipeline.

The `multiqc_config.yml` in `assets/` sets:
- Custom report title and subtitle
- Which modules to include/exclude
- Sample name cleanup regexes (removes `.fastq.gz`, `.bam` etc. from names)

---

## Assets: samplesheet.csv and metadata.csv

### samplesheet.csv

```
sample,condition,fastq1,fastq2
SRR1039508,control,data/SRR1039508_1.fastq.gz,data/SRR1039508_2.fastq.gz
...
```

Parsed in `main.nf`. The `sample` field becomes:
- The `meta.id` in all Nextflow channels (used in log tags, output file names)
- The column name in `gene_counts.txt` (after `.bam` stripping)
- The `rownames` of the DESeq2 metadata

### metadata.csv

```
sample,condition,cell_line,batch,sex
SRR1039508,control,N61311,1,male
...
```

Used as the `colData` of the DESeqDataSet. Each column becomes a variable available in the DESeq2 design formula. Currently only `condition` is used in the `~ condition` design, but `batch`, `cell_line`, and `sex` are present for researchers who want to model covariates (e.g., `~ batch + condition`).

The 4 cell lines (N61311, N052611, N080611, N061011) are 4 different human donors — biological replicates. Each appears once as control and once as treated, making this a **paired design**. This dramatically increases statistical power.

---

## How to interpret each output file

### Is my sample good quality?

Open `results/multiqc/multiqc_report.html`. Look at:

1. **FastQC — Per base sequence quality**: Green = good. Boxes (median) should be in the green zone (Q > 28) for all cycles.
2. **fastp — After filtering**: Read length distribution should still peak near 150 bp (not heavily trimmed).
3. **STAR — % Uniquely Mapped**: Should be > 88%. Red if < 80%.
4. **featureCounts — Assigned**: Should be > 60%. Low numbers suggest a GTF mismatch.
5. **samtools — % Properly Paired**: Should be > 90%.

### Are there real DEGs?

Open `results/deseq2/plots/Volcano_plot.png`:
- If you see a clear symmetrical split with many red (up) and blue (down) dots, the treatment has a strong transcriptional effect.
- If the plot is flat (all grey, nothing above the p-value line), either the treatment had no effect, the sample QC was poor, or the thresholds are too strict.

Open `results/deseq2/plots/PCA_plot.png`:
- If control and treated samples form separate clusters, you have real biological signal.
- If samples cluster by cell line instead, consider adding `cell_line` or `batch` to the DESeq2 design.

### What are the most interesting genes?

Open `results/deseq2/DEG_significant.csv`. Sort by `padj` ascending (most significant first). The top hits for this dataset will be:

| Gene | Direction | Biology |
|------|-----------|---------|
| FKBP5 | Up | GR co-chaperone; classic dexamethasone target |
| CRISPLD2 | Up | Anti-inflammatory secreted protein |
| DUSP1 | Up | Dual-specificity phosphatase; suppresses MAPK |
| TSC22D3 (GILZ) | Up | Transcription factor; mediates GR anti-inflammatory effects |
| IL8 (CXCL8) | Down | Pro-inflammatory cytokine |
| CXCL1 | Down | Neutrophil chemokine |
| CCL2 | Down | Monocyte chemokine |

### What pathways are changed?

Open `results/pathway_analysis/reactome_ora/reactome_ORA_results.csv`. Sort by `p.adjust` ascending. For GSE52778, the top pathways will typically be:

**Activated (up-regulated DEGs):**
- Glucocorticoid receptor signalling
- Gene expression (Transcription)
- mRNA processing

**Suppressed (down-regulated DEGs):**
- Interleukin-1 signalling
- TNF signalling
- MAPK cascade
- Innate immune system
- Cytokine signalling in immune system

For the GSEA view, open `results/pathway_analysis/reactome_gsea/reactome_GSEA_results.csv`. NES > 0 = activated, NES < 0 = suppressed.

### Can I reuse the results in R?

Yes. Load the saved R objects:

```r
library(DESeq2)

# Load the full DESeqDataSet (model, counts, colData)
dds <- readRDS("results/deseq2/rds/dds.rds")

# Load results
res <- readRDS("results/deseq2/rds/res.rds")

# Load VST-transformed data
vsd <- readRDS("results/deseq2/rds/vsd.rds")

# Example: custom plot for a specific gene
plotCounts(dds, gene="ENSG00000096060", intgroup="condition")  # FKBP5
```

These `.rds` files preserve every aspect of the analysis — you can run any DESeq2 function on them.

---

## Glossary

| Term | Definition |
|------|-----------|
| **BAM** | Binary Alignment Map — compressed binary format storing each read's alignment position |
| **BAI** | BAM index — allows random access to a BAM by genomic coordinate |
| **baseMean** | Average normalised read count across all samples for a gene |
| **BH correction** | Benjamini-Hochberg False Discovery Rate correction for multiple testing |
| **DEG** | Differentially Expressed Gene — a gene whose expression changes significantly between conditions |
| **Dispersion** | DESeq2's measure of gene-wise variability beyond what the mean would predict |
| **FASTQ** | Text format for raw sequencing reads: sequence + quality scores |
| **FDR** | False Discovery Rate — expected proportion of false positives among significant results |
| **GSEA** | Gene Set Enrichment Analysis — tests whether a gene set is enriched at one end of a ranked list |
| **GTF** | Gene Transfer Format — annotation file listing gene/exon coordinates |
| **GR** | Glucocorticoid Receptor |
| **LFC** | Log2 Fold Change — log2(treated/control); +1 = 2-fold up, -1 = 2-fold down |
| **NES** | Normalized Enrichment Score — the GSEA effect size statistic |
| **ORA** | Over-Representation Analysis — tests if a pathway is over-represented in a DEG list |
| **padj** | Adjusted p-value (BH-corrected) |
| **PCA** | Principal Component Analysis — dimensionality reduction method |
| **Phred score** | Quality score: Q30 = 1 in 1000 chance of incorrect base call |
| **VST** | Variance Stabilising Transformation — makes counts comparable across expression levels |
