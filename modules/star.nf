// modules/star.nf
// STAR alignment with piRNA-specific parameters
// Faithfully translated from your original loop:
//   --outFilterMismatchNmax 0
//   --outFilterMultimapNmax 1000000
//   --outFilterMatchNmin 22
//   --outFilterMultimapScoreRange 1
//   --alignEndsType Local
//   --outSAMunmapped Within
//   --outReadsUnmapped None

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.out.bam"), emit: bam
    path("${sample_id}_Log.final.out"),                          emit: log_final
    path("${sample_id}_Log.out"),                                emit: log_out

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode alignReads \\
        --genomeDir ${params.genome_dir} \\
        --readFilesIn ${trimmed_reads} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes Standard \\
        --outSAMunmapped Within \\
        --outReadsUnmapped None \\
        --alignEndsType Local \\
        --outFilterMismatchNmax ${params.star_mismatch} \\
        --outFilterMultimapScoreRange 1 \\
        --outFilterMatchNmin ${params.star_min_match} \\
        --outFilterMultimapNmax ${params.star_multimap} \\
        --sjdbGTFfeatureExon ${params.gtf} \\
        --outFileNamePrefix ${sample_id}_
    """

    // STAR names output as <prefix>Aligned.out.bam by default
    // We use sample_id as prefix so outputs are clearly labelled
}
