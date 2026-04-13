// modules/fastp.nf
// Fastp: generalized adapter trimming for multi-dataset piRNA analysis
//
// Key design: NO hardcoded adapter sequence
//   - fastp auto-detects adapters from all Illumina kits
//     (TruSeq, DpnII, NEBNext, SMARTer, etc.)
//   - --trim_poly_g handles NextSeq/NovaSeq 2-color artifacts (no-op on HiSeq)
//   - --trim_poly_x handles polyA tailing from library prep (no-op if absent)
//   - Length filter standardizes output across all read lengths (50-150bp input)

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp",        mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.outdir}/fastp/reports", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fastp.fastq.gz"), emit: trimmed
    path("${sample_id}.fastp.html"),                           emit: html
    path("${sample_id}.fastp.json"),                           emit: json

    script:
    """
    fastp \\
        -i ${reads} \\
        -o ${sample_id}.fastp.fastq.gz \\
        --trim_poly_g \\
        --poly_g_min_len ${params.poly_g_min_len} \\
        --trim_poly_x \\
        --poly_x_min_len ${params.poly_x_min_len} \\
        --length_required ${params.min_len} \\
        --length_limit ${params.max_len} \\
        --qualified_quality_phred ${params.quality_phred} \\
        --low_complexity_filter \\
        --complexity_threshold 30 \\
        --thread ${task.cpus} \\
        --html ${sample_id}.fastp.html \\
        --json ${sample_id}.fastp.json
    """
}
