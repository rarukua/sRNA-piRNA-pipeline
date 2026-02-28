// modules/fastqc.nf
// Run FastQC on raw input reads

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}
    """
}
