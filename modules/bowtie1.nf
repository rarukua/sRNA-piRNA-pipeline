// modules/bowtie1.nf
// Bowtie1 alignment — optimized for piRNA short reads
//
// -v 1          : allow 1 mismatch (captures SNV-containing isoforms)
// -a            : report ALL alignments in best stratum (no cap)
//                 piRNAs from TEs can map to >100 loci
// --best        : report best alignments only
// --strata      : only report alignments in best stratum
// -S            : output SAM format
// --chunkmbs    : memory per thread for multi-mapper resolution

process BOWTIE1_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/bowtie1", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    path("${sample_id}.bowtie1.log"),                          emit: log
    path("${sample_id}.flagstat"),                             emit: flagstat

    script:
    """
    bowtie \\
        -v ${params.bowtie1_mismatch} \\
        -a \\
        --best \\
        --strata \\
        --sam \\
        -p ${task.cpus} \\
        --chunkmbs 512 \\
        -x ${params.bowtie1_index} \\
        <(zcat ${trimmed_reads}) \\
        2> ${sample_id}.bowtie1.log \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat
    """
}
