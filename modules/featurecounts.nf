// modules/featurecounts.nf
// featureCounts quantification — piRNA counting with multi-mapper + overlap handling
//
// Pipeline: filter BAM to piRNA-sized reads (≥24nt) → featureCounts
//
// Key flags:
//   -t exon        : count against 'exon' lines in GTF (piRNA countable regions)
//   -M --fraction  : multi-mappers counted fractionally (1/N per locus)
//   -O --fraction  : overlapping features counted fractionally (1/N per feature)
//                    Essential for dense piRNA clusters where isoforms overlap
//   -Q 0           : no MAPQ filter (bowtie1 multi-mappers have MAPQ 0)
//   -s 1           : forward-stranded (piRNAs are strand-specific)
//   --minOverlap   : minimum read-feature overlap (15nt balances specificity vs isoform capture)
//   --largestOverlap : tiebreak by largest overlap

process FEATURECOUNTS {
    tag "$sample_id"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_counts.txt"),         emit: counts
    path("${sample_id}_counts.txt.summary"),                        emit: summary

    script:
    """
    # Filter BAM to piRNA-sized reads (≥24nt) to remove rRNA/tRNA junk
    samtools view -h ${bam} | \\
        awk 'BEGIN{OFS="\\t"} /^@/ || length(\$10)>=24' | \\
        samtools view -bS -o ${sample_id}.pirna_sized.bam -
    samtools index -@ ${task.cpus} ${sample_id}.pirna_sized.bam

    # Count with fractional multi-mapper + overlap handling
    featureCounts \\
        -t ${params.fc_feature_type} \\
        -F GTF \\
        -a ${params.gtf} \\
        -M \\
        -O \\
        --fraction \\
        --minOverlap ${params.min_overlap} \\
        --largestOverlap \\
        -Q 0 \\
        -s ${params.fc_strandedness} \\
        -T ${task.cpus} \\
        -o ${sample_id}_counts.txt \\
        ${sample_id}.pirna_sized.bam

    # Clean up intermediate
    rm -f ${sample_id}.pirna_sized.bam ${sample_id}.pirna_sized.bam.bai
    """
}
