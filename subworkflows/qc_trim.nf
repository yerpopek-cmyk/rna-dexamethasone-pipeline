/*
================================================================================
    QC_TRIM - Subworkflow: Quality control and adapter trimming
================================================================================
*/

include { FASTQC  } from '../modules/fastqc'
include { FASTP   } from '../modules/fastp'

workflow QC_TRIM {

    take:
    reads  // channel: [meta, r1, r2]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Pre-trimming FastQC
    //
    FASTQC(reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Adapter trimming with fastp
    //
    FASTP(reads)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    emit:
    trimmed_reads  = FASTP.out.trimmed   // channel: [meta, r1_trimmed, r2_trimmed]
    fastqc_zip     = FASTQC.out.zip      // channel: [meta, *.zip]
    fastqc_html    = FASTQC.out.html     // channel: [meta, *.html]
    trim_log       = FASTP.out.log       // channel: [*.log]
    trim_html      = FASTP.out.html      // channel: [*.html]
    trim_json      = FASTP.out.json      // channel: [*.json]
    versions       = ch_versions         // channel: [versions.yml]
}
