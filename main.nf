#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
================================================================================
    rna-dexamethasone: RNA-seq Analysis Pipeline
    Dexamethasone effects on airway smooth muscle cells (GSE52778)
================================================================================
    Author      : Bioinformatics Pipeline
    Version     : 1.0.0
    Nextflow    : >=23.04.0
    Description : Complete RNA-seq workflow from FASTQ to pathway enrichment
================================================================================
*/

// Import subworkflows
include { QC_TRIM          } from './subworkflows/qc_trim'
include { ALIGN_QUANTIFY   } from './subworkflows/align_quantify'
include { DIFF_EXPRESSION  } from './subworkflows/differential_expression'

// Import standalone modules
include { MULTIQC          } from './modules/multiqc'

// Validate parameters
def validateParams() {
    if (!params.samplesheet) {
        error "ERROR: --samplesheet is required. Please provide a CSV samplesheet."
    }
    if (!params.star_index && !params.genome) {
        error "ERROR: Either --star_index or --genome must be specified."
    }
    if (!params.gtf) {
        error "ERROR: --gtf (genome annotation file) is required."
    }
}

// Print pipeline header
def printHeader() {
    log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║              RNA-seq Dexamethasone Analysis Pipeline  v${workflow.manifest.version}              ║
║              Airway Smooth Muscle Cells (GSE52778)                         ║
╚══════════════════════════════════════════════════════════════════════════════╝

Pipeline parameters:
──────────────────────────────────────────────────────────
  Samplesheet    : ${params.samplesheet}
  Metadata       : ${params.metadata}
  Genome         : ${params.genome}
  STAR index     : ${params.star_index}
  GTF            : ${params.gtf}
  Output dir     : ${params.outdir}
  p-value cutoff : ${params.pvalue_cutoff}
  LFC threshold  : ${params.lfc_threshold}
  Organism       : ${params.organism}
──────────────────────────────────────────────────────────
  Run name       : ${workflow.runName}
  Profile        : ${workflow.profile}
  Container      : ${workflow.container ?: 'N/A'}
──────────────────────────────────────────────────────────
    """.stripIndent()
}

// ============================================================
//   MAIN WORKFLOW
// ============================================================

workflow {
    main:
    validateParams()
    printHeader()

    // ─── Input channel from samplesheet ───────────────────────
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id        : row.sample,
                condition : row.condition,
                single_end: row.fastq2 ? false : true
            ]
            def r1 = file(row.fastq1, checkIfExists: true)
            def r2 = row.fastq2 ? file(row.fastq2, checkIfExists: true) : []
            return [meta, r1, r2]
        }
        .set { ch_reads }

    // ─── Metadata for R analysis ──────────────────────────────
    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)

    // ─── SUBWORKFLOW: Quality control and trimming ─────────────
    QC_TRIM(ch_reads)

    // ─── SUBWORKFLOW: Alignment and quantification ─────────────
    ALIGN_QUANTIFY(QC_TRIM.out.trimmed_reads)

    // ─── SUBWORKFLOW: Differential expression + pathway analysis
    DIFF_EXPRESSION(
        ALIGN_QUANTIFY.out.counts,
        ch_metadata
    )

    // ─── MODULE: MultiQC aggregation ──────────────────────────
    ch_multiqc_files = Channel.empty()
        .mix(QC_TRIM.out.fastqc_zip.map{ meta, files -> files }.collect())
        .mix(QC_TRIM.out.trim_log.collect())
        .mix(ALIGN_QUANTIFY.out.star_log.map{ meta, file -> file }.collect())
        .mix(ALIGN_QUANTIFY.out.featurecounts_summary.collect())
        .mix(ALIGN_QUANTIFY.out.flagstat.map{ meta, file -> file }.collect())
        .mix(ALIGN_QUANTIFY.out.bam_stats.map{ meta, file -> file }.collect())
        .mix(ALIGN_QUANTIFY.out.idxstats.map{ meta, file -> file }.collect())
        .collect()

    MULTIQC(ch_multiqc_files)

}

// ============================================================
//   WORKFLOW COMPLETION
// ============================================================
