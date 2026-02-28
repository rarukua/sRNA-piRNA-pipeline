// modules/bowtie1.nf
// Bowtie1 alignment — optimized for piRNA short reads
// -v 0          : 0 mismatches (matches your STAR outFilterMismatchNmax 0)
// --best        : report best alignments only
// --strata      : only report alignments in best stratum
// -a            : report all alignments in best stratum (like STAR's 1M multimappers)
// -S            : output SAM format

process BOWTIE1_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/bowtie1", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path("${sample_id}.bowtie1.log"),               emit: log

    script:
    """
    # Decompress first since bowtie1 doesn't support zcat pipe natively
    gunzip -c ${trimmed_reads} > ${sample_id}.fastp.fastq

    bowtie \
        -p ${task.cpus} \
        -v 0 \
        --best \
        --strata \
        -a \
        -x ${params.bowtie1_index} \
        ${sample_id}.fastp.fastq \
        -S ${sample_id}.sam \
        2> ${sample_id}.bowtie1.log

    # Convert SAM to BAM and clean up
    samtools view -bS ${sample_id}.sam | samtools sort -o ${sample_id}.bam
    rm ${sample_id}.sam ${sample_id}.fastp.fastq
    """
}
