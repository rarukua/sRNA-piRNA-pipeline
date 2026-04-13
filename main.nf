#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """
    ╔═══════════════════════════════════════════════╗
    ║         sRNA / piRNA Analysis Pipeline        ║
    ║              v2.0 — Generalized               ║
    ╚═══════════════════════════════════════════════╝
    samplesheet   : ${params.samplesheet}
    aligner       : ${params.aligner}
    bowtie1 index : ${params.bowtie1_index}
    annotation    : ${params.gtf}
    min_len       : ${params.min_len}
    max_len       : ${params.max_len}
    min_overlap   : ${params.min_overlap}
    mode          : ${params.mode}
    outdir        : ${params.outdir}
    """.stripIndent()

include { FASTQC                   } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { FASTP                    } from './modules/fastp'
include { STAR_ALIGN               } from './modules/star'
include { BOWTIE2_ALIGN            } from './modules/bowtie2'
include { BOWTIE1_ALIGN            } from './modules/bowtie1'
include { FEATURECOUNTS            } from './modules/featurecounts'
include { MERGE_COUNTS             } from './modules/merge_counts'
include { MULTIQC                  } from './modules/multiqc'

workflow {

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq)) }
        .set { reads_ch }

    // ── QC on raw reads ──
    FASTQC(reads_ch)

    // ── Trimming: auto-detect adapters + polyG/polyX + size selection ──
    FASTP(reads_ch)
    FASTQC_TRIMMED(FASTP.out.trimmed)

    if (params.mode == 'qc_only') {

        log.info "QC-only mode — check MultiQC report before proceeding to alignment"

        qc_only_ch = FASTQC.out.zip.map { it[1] }
            .mix( FASTQC_TRIMMED.out.zip.map { it[1] } )
            .mix( FASTP.out.json )
            .collect()

        MULTIQC(qc_only_ch)

    } else {

        // ── Alignment ──
        if (params.aligner == 'star') {
            STAR_ALIGN(FASTP.out.trimmed)
            bam_ch = STAR_ALIGN.out.bam
        } else if (params.aligner == 'bowtie2') {
            BOWTIE2_ALIGN(FASTP.out.trimmed)
            bam_ch = BOWTIE2_ALIGN.out.bam
        } else if (params.aligner == 'bowtie1') {
            BOWTIE1_ALIGN(FASTP.out.trimmed)
            bam_ch = BOWTIE1_ALIGN.out.bam
        } else {
            error "Invalid aligner '${params.aligner}'. Choose 'star', 'bowtie2', or 'bowtie1'."
        }

        // ── Quantification: filter ≥24nt + featureCounts ──
        FEATURECOUNTS(bam_ch)

        // ── Merge count tables into matrices ──
        MERGE_COUNTS(
            FEATURECOUNTS.out.counts.map { it[1] }.collect()
        )

        // ── MultiQC: aggregate all QC ──
        qc_files_ch = FASTQC.out.zip.map { it[1] }
            .mix( FASTQC_TRIMMED.out.zip.map { it[1] } )
            .mix( FASTP.out.json )
            .mix( FEATURECOUNTS.out.summary )
            .collect()

        MULTIQC(qc_files_ch)
    }
}
