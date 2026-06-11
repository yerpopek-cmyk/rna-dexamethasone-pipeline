/*
================================================================================
    ALIGN_QUANTIFY - Subworkflow: Alignment and read quantification
================================================================================
*/

include { STAR_ALIGN    } from '../modules/star_align'
include { BAM_QC        } from '../modules/bam_qc'
include { FEATURECOUNTS } from '../modules/featurecounts'

workflow ALIGN_QUANTIFY {

    take:
    reads  // channel: [meta, r1_trimmed, r2_trimmed]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Align reads to reference genome with STAR
    //
    STAR_ALIGN(reads)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    // Join BAM and BAI
    ch_bam_bai = STAR_ALIGN.out.bam
        .join(STAR_ALIGN.out.bai)

    //
    // MODULE: BAM quality control (consolidated)
    //
    BAM_QC(ch_bam_bai)
    ch_versions = ch_versions.mix(BAM_QC.out.versions.first())

    //
    // MODULE: Count reads per gene
    //
    // Collect all BAM files for featureCounts (single call for the count matrix)
    // Convert to absolute path strings to bypass file staging/copying in Nextflow
    ch_bam_all = STAR_ALIGN.out.bam
        .map { meta, bam -> bam.toAbsolutePath().toString() }
        .collect()

    FEATURECOUNTS(ch_bam_all, params.gtf)
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)

    emit:
    bam                   = STAR_ALIGN.out.bam          // channel: [meta, *.bam]
    bai                   = STAR_ALIGN.out.bai          // channel: [meta, *.bai]
    star_log              = STAR_ALIGN.out.log_final    // channel: [meta, *Log.final.out]
    flagstat              = BAM_QC.out.flagstat         // channel: [meta, *.flagstat]
    bam_stats             = BAM_QC.out.stats            // channel: [meta, *.stats]
    idxstats              = BAM_QC.out.idxstats         // channel: [meta, *.idxstats]
    counts                = FEATURECOUNTS.out.counts    // path: gene_counts.txt
    featurecounts_summary = FEATURECOUNTS.out.summary  // path: *.summary
    versions              = ch_versions
}
