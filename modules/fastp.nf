// modules/fastp.nf
// Fastp: adapter trimming + piRNA size selection
// Your original: fastp -i $file -o ... --max_len1 24
// We preserve all your original params and expose them via nextflow.config

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
        --max_len1 ${params.max_len} \\
        --thread ${task.cpus} \\
        --html ${sample_id}.fastp.html \\
        --json ${sample_id}.fastp.json
    """
}

// NOTE on your original params:
// --max_len1 24     → keeps only reads ≤ 24 nt (piRNA size selection)
// No --min_len set  → add --min_len 18 in nextflow.config if you want a lower bound
// No adapter set    → fastp auto-detects adapters (fine for most cases)
