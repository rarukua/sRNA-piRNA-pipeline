// modules/bowtie2.nf
// Bowtie2 alignment — piRNA-appropriate parameters
// piRNA-equivalent settings:
//   --no-unal         → suppress unaligned reads in output (like STAR's outReadsUnmapped None)
//   --local           → local alignment (like STAR's alignEndsType Local)
//   --score-min L,22,0 → minimum alignment score = 22 (like STAR's outFilterMatchNmin 22)
//   -N 0              → no mismatches in seed (like STAR's outFilterMismatchNmax 0)
//   -k 1000000        → report up to 1M alignments per read (like STAR's outFilterMultimapNmax)

process BOWTIE2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path("${sample_id}.bowtie2.log"),               emit: log

    script:
    """
    bowtie2 \\
        -p ${task.cpus} \\
        -x ${params.bowtie2_index} \\
        -U ${trimmed_reads} \\
        --local \\
        -N 0 \\
        --score-min L,${params.min_overlap},0 \\
        -k 1000000 \\
        --no-unal \\
        2> ${sample_id}.bowtie2.log \\
    | samtools view -bS -o ${sample_id}.bam
    """

    // NOTE: Bowtie2 doesn't natively output BAM — we pipe through samtools view
    // samtools is included in the Docker image
}
