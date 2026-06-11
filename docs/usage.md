# Usage Guide

## Prerequisites

- Nextflow >= 23.04.0
- Docker or Conda
- 64 GB RAM (for STAR alignment)
- 8+ CPU cores recommended

## Quick start

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Clone the pipeline
```bash
git clone https://github.com/yerpopek-cmyk/rna-dexamethasone-pipeline.git
cd rna-dexamethasone
```

### 3. Download reference genome (GRCh38)
```bash
# Download FASTA and GTF from GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

# Build STAR index (~30 min, 32 GB RAM)
mkdir -p refs/star_index
STAR --runMode genomeGenerate \
     --genomeDir refs/star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v46.annotation.gtf \
     --runThreadN 16 \
     --genomeSAindexNbases 14
```

### 4. Prepare FASTQ data (GSE52778)
```bash
# Download SRA data using SRA toolkit
for SRR in SRR1039508 SRR1039509 SRR1039512 SRR1039513 \
            SRR1039516 SRR1039517 SRR1039520 SRR1039521; do
    fasterq-dump --split-files --gzip --outdir data/ $SRR
done
```

### 5. Edit samplesheet
Edit `assets/samplesheet.csv` to point to your FASTQ files:
```csv
sample,condition,fastq1,fastq2
SRR1039508,control,data/SRR1039508_1.fastq.gz,data/SRR1039508_2.fastq.gz
```

### 6. Run the pipeline

**With Docker:**
```bash
nextflow run main.nf \
    -profile docker \
    --star_index refs/star_index \
    --gtf gencode.v46.annotation.gtf \
    --outdir results
```

**With Conda:**
```bash
nextflow run main.nf \
    -profile conda \
    --star_index refs/star_index \
    --gtf gencode.v46.annotation.gtf \
    --outdir results
```

**On SLURM HPC:**
```bash
nextflow run main.nf \
    -profile slurm,docker \
    --star_index refs/star_index \
    --gtf gencode.v46.annotation.gtf \
    --outdir results \
    -resume
```

## Using the airway R package (skip alignment)

If you want to start directly from the pre-computed count matrix in R,
run the analysis scripts manually:

```bash
Rscript bin/deseq2.R \
    --counts  results/featurecounts/gene_counts.txt \
    --metadata assets/metadata.csv \
    --outdir   results/deseq2

Rscript bin/reactome.R \
    --deg     results/deseq2/DEG_significant.csv \
    --outdir  results/pathway_analysis/reactome_ora

Rscript bin/gsea.R \
    --deg     results/deseq2/DEG_results.csv \
    --outdir  results/pathway_analysis/reactome_gsea
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `assets/samplesheet.csv` | Input samplesheet |
| `--metadata` | `assets/metadata.csv` | Sample metadata for R |
| `--star_index` | `null` | Path to STAR genome index |
| `--gtf` | `null` | Genome annotation GTF |
| `--outdir` | `results` | Output directory |
| `--pvalue_cutoff` | `0.05` | Adjusted p-value threshold |
| `--lfc_threshold` | `1.0` | \|log2FC\| threshold for DEGs |
| `--organism` | `human` | Organism for ReactomePA |
| `--deseq2_min_count` | `10` | Minimum total count filter |
| `--star_threads` | `8` | Threads for STAR alignment |

## Output structure

```
results/
├── fastqc/              # Per-sample FastQC HTML/ZIP
├── fastp/               # Trimming reports (HTML, JSON)
├── star/                # BAM files + STAR logs
├── bam_qc/              # samtools flagstat/stats
├── featurecounts/       # gene_counts.txt + summary
├── deseq2/
│   ├── DEG_results.csv
│   ├── DEG_significant.csv
│   ├── normalized_counts.csv
│   ├── vst_counts.csv
│   ├── plots/           # PCA, Volcano, Heatmap, MA, ...
│   └── rds/             # DESeq2 R objects
├── pathway_analysis/
│   ├── reactome_ora/    # ORA results + plots
│   └── reactome_gsea/   # GSEA results + plots
├── multiqc/             # MultiQC report
├── report/              # Final HTML report
└── pipeline_info/       # Execution logs, timeline, DAG
```
