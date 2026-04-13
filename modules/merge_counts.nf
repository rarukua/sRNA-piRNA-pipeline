// modules/merge_counts.nf
// Merge individual featureCounts outputs into combined matrices
// Produces: raw count matrix, RPM-normalized matrix, library size summary

process MERGE_COUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(count_files)  // collected list of all *_counts.txt files

    output:
    path("piRNA_counts_matrix.tsv"),    emit: raw_matrix
    path("piRNA_rpm_matrix.tsv"),       emit: rpm_matrix
    path("library_sizes.tsv"),          emit: lib_sizes

    script:
    """
    #!/usr/bin/env python3
    import os, glob

    count_files = sorted(glob.glob("*_counts.txt"))
    if not count_files:
        raise RuntimeError("No count files found!")

    gene_ids = []
    gene_info = {}
    sample_counts = {}

    for fi, cf in enumerate(count_files):
        sample = os.path.basename(cf).replace("_counts.txt", "")
        counts = {}

        with open(cf) as f:
            for line in f:
                if line.startswith('#') or line.startswith('Geneid'):
                    continue
                parts = line.strip().split('\\t')
                gene_id = parts[0]
                if fi == 0:
                    gene_ids.append(gene_id)
                    gene_info[gene_id] = parts[1:6]
                count = float(parts[-1])
                counts[gene_id] = count

        sample_counts[sample] = counts

    samples = sorted(sample_counts.keys())

    # Raw counts matrix
    with open("piRNA_counts_matrix.tsv", 'w') as f:
        f.write("gene_id\\t" + "\\t".join(samples) + "\\n")
        for gid in gene_ids:
            row = [gid] + [f"{sample_counts[s].get(gid, 0):.2f}" for s in samples]
            f.write("\\t".join(row) + "\\n")

    # Library sizes
    lib_sizes = {s: sum(sample_counts[s].values()) for s in samples}

    with open("library_sizes.tsv", 'w') as f:
        f.write("sample\\ttotal_assigned_reads\\n")
        for s in samples:
            f.write(f"{s}\\t{lib_sizes[s]:.2f}\\n")

    # RPM-normalized matrix
    with open("piRNA_rpm_matrix.tsv", 'w') as f:
        f.write("gene_id\\t" + "\\t".join(samples) + "\\n")
        for gid in gene_ids:
            row = [gid]
            for s in samples:
                raw = sample_counts[s].get(gid, 0)
                rpm = (raw / lib_sizes[s]) * 1_000_000 if lib_sizes[s] > 0 else 0
                row.append(f"{rpm:.4f}")
            f.write("\\t".join(row) + "\\n")

    print(f"Merged {len(count_files)} samples, {len(gene_ids)} features")
    """
}
