#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """
    ╔═══════════════════════════════════════════════╗
    ║         sRNA / piRNA Analysis Pipeline        ║
    ╚═══════════════════════════════════════════════╝
    samplesheet   : ${params.samplesheet}
    genome_dir    : ${params.genome_dir}
    annotation    : ${params.gtf}
    aligner       : ${params.aligner}
    max_len       : ${params.max_len}
    min_overlap   : ${params.min_overlap}
    mode          : ${params.mode}
    outdir        : ${params.outdir}
    """.stripIndent()

include { FASTQC              } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { FASTP               } from './modules/fastp'
include { STAR_ALIGN          } from './modules/star'
include { BOWTIE2_ALIGN       } from './modules/bowtie2'
include { FEATURECOUNTS       } from './modules/featurecounts'
include { MULTIQC             } from './modules/multiqc'

workflow {

    // --- 1. Read samplesheet ---
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq)) }
        .set { reads_ch }

    // --- 2. FastQC on raw reads ---
    FASTQC(reads_ch)

    // --- 3. Fastp trimming + size selection ---
    FASTP(reads_ch)

    // --- 4. FastQC on trimmed reads (so you can compare before/after) ---
    FASTQC_TRIMMED(FASTP.out.trimmed)

    // ================================================================
    // MODE: qc_only
    // Stops here and runs MultiQC on raw + trimmed FastQC results.
    // Use this FIRST on a new dataset to validate your trim parameters
    // before committing to alignment.
    //
    // Run with: --mode qc_only
    // ================================================================
    if (params.mode == 'qc_only') {

        log.info "Running in QC-only mode — stopping after FastQC + fastp"
        log.info "Check results/multiqc/multiqc_report.html to validate trim params"
        log.info "Then rerun with --mode full to proceed to alignment"

        qc_only_ch = FASTQC.out.zip.map { it[1] }
            .mix( FASTQC_TRIMMED.out.zip.map { it[1] } )
            .mix( FASTP.out.json )
            .collect()

        MULTIQC(qc_only_ch)

    // ================================================================
    // MODE: full (default)
    // Runs the complete pipeline: trim → align → quantify → MultiQC
    // Only use after validating QC params with qc_only mode first.
    // ================================================================
    } else {

        if (params.aligner == 'star') {
            STAR_ALIGN(FASTP.out.trimmed)
            bam_ch = STAR_ALIGN.out.bam
        } else if (params.aligner == 'bowtie2') {
            BOWTIE2_ALIGN(FASTP.out.trimmed)
            bam_ch = BOWTIE2_ALIGN.out.bam
        } else {
            error "Invalid aligner '${params.aligner}'. Choose 'star' or 'bowtie2'."
        }

        FEATURECOUNTS(bam_ch)

        // Final MultiQC includes everything: raw QC, trimmed QC, alignment, counts
        qc_files_ch = FASTQC.out.zip.map { it[1] }
            .mix( FASTQC_TRIMMED.out.zip.map { it[1] } )
            .mix( FASTP.out.json )
            .mix( FEATURECOUNTS.out.summary )
            .collect()

        MULTIQC(qc_files_ch)
    }
}
