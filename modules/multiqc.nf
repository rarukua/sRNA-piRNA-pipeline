// modules/multiqc.nf
// MultiQC: aggregate QC from FastQC, fastp, STAR/Bowtie2, featureCounts

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(qc_files)  // collected list of all QC output files

    output:
    path("*_multiqc_report.html"), emit: report

    script:
    """
    multiqc \\
        --title "${params.multiqc_title}" \\
        --force \\
        --outdir . \\
        .
    """
}
