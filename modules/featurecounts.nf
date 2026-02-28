// modules/featurecounts.nf
// featureCounts quantification — piRNA cluster counting
// From your original:
//   featureCounts -t piRbase -F GTF
//                -a piRbase_cluster_GTF.gtf
//                -M -O --minOverlap 22 -Q 0

process FEATURECOUNTS {
    tag "$sample_id"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_cluster_counts.txt"),         emit: counts
    path("${sample_id}_cluster_counts.txt.summary"),                        emit: summary

    script:
    """
    featureCounts \\
        -t ${params.fc_feature_type} \\
        -F GTF \\
        -a ${params.gtf} \\
        -M \\
        -O \\
        --minOverlap ${params.min_overlap} \\
        -Q 0 \\
        -T ${task.cpus} \\
        -o ${sample_id}_cluster_counts.txt \\
        ${bam}
    """

    // Key piRNA flags explained:
    // -t piRbase       → feature type matching your GTF 3rd column
    // -M               → count multi-mapping reads (critical for piRNAs — they map to repeats)
    // -O               → allow reads to be assigned to overlapping features
    // --minOverlap 22  → minimum overlap between read and feature = 22 nt
    // -Q 0             → no minimum mapping quality filter (needed for multimappers)
}
