# sRNA / piRNA Analysis Pipeline

A reproducible Nextflow pipeline for small RNA-seq analysis with a focus on **piRNA biology**. Supports STAR, Bowtie2, and Bowtie1 alignment, fastp-based size selection, and piRNA cluster quantification with featureCounts.

---

## Pipeline Overview

```
Raw FASTQ
    │
    ▼
 FastQC (raw)        ← QC on raw reads
    │
    ▼
 Fastp               ← Adapter trimming + size selection (≤24 nt for piRNA)
    │
    ▼
 FastQC (trimmed)    ← QC after trimming (compare before/after in MultiQC)
    │
    ▼
 STAR                ← Splice-aware, multimapper-aware (fastest, recommended)
   or
 Bowtie1             ← Short-read optimized, piRNA community standard
   or
 Bowtie2             ← General-purpose aligner (slower for piRNA)
    │
    ▼
 featureCounts       ← piRNA cluster quantification against piRBase GTF
    │
    ▼
 MultiQC             ← Aggregated QC report (raw + trimmed + alignment + counts)
```

---

## Two-Stage Recommended Workflow

This pipeline supports two modes designed to match real-world analysis practice:

### Stage 1 — QC only (run first on any new dataset)
```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --mode qc_only \
    --outdir results_qc
```
Opens `results_qc/multiqc/` — shows raw vs trimmed FastQC side by side for all samples.
Verify adapter removal worked and read length peaks at ~24 nt before proceeding.

### Stage 2 — Full pipeline (after QC validation)
```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --mode full \
    --aligner star \          # or bowtie1 / bowtie2
    --genome_dir /path/to/STAR_index \
    --gtf /path/to/piRbase_cluster.gtf \
    --outdir results_full \
    -resume                   # skips already-completed steps
```

---

## Key piRNA-Specific Parameters

| Step | Parameter | Value | Rationale |
|------|-----------|-------|-----------|
| fastp | `--max_len1` | 24 | Retain reads ≤24 nt (piRNA size range) |
| STAR | `--outFilterMismatchNmax` | 0 | Exact matches only |
| STAR | `--outFilterMultimapNmax` | 1,000,000 | Allow multimappers (piRNAs map to repeats) |
| STAR | `--outFilterMatchNmin` | 22 | Minimum matched bases |
| STAR | `--alignEndsType` | Local | Soft-clipping allowed |
| Bowtie1 | `-v` | 0 | 0 mismatches |
| Bowtie1 | `--best --strata -a` | enabled | All best-stratum alignments |
| Bowtie2 | `-N` | 0 | No mismatches in seed |
| Bowtie2 | `-k` | 1,000,000 | Report up to 1M alignments |
| featureCounts | `--minOverlap` | 22 | Minimum read-feature overlap |
| featureCounts | `-M -O` | enabled | Count multimappers + overlapping features |
| featureCounts | `-Q` | 0 | No MAPQ filter (needed for multimappers) |
| featureCounts | `-t` | piRbase | Feature type matching piRBase GTF |

---

## Aligner Comparison

| Aligner | Speed | Best for | Notes |
|---------|-------|----------|-------|
| STAR | ⚡ Fastest | General piRNA-seq | Splice-aware, handles multimappers efficiently |
| Bowtie1 | 🐢 Moderate | piRNA (community standard) | Designed for short reads ≤50nt, used by proTRAC/piPipes |
| Bowtie2 | 🐢 Slowest | General sRNA | `-k 1000000` makes it slow for repetitive piRNA loci |

**Recommendation:** Use STAR for speed, Bowtie1 for compatibility with piRNA-specific downstream tools.

---

## Installation

### 1. Install Java 17+ and Nextflow
```bash
# Install SDKMAN (no sudo needed)
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"

# Install Java 17
sdk install java 17.0.10-tem

# Install Nextflow
mkdir -p ~/bin
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

nextflow -version
```

### 2. Install required tools (via conda)
```bash
conda install -c bioconda -c conda-forge \
    fastqc fastp star bowtie bowtie2 samtools subread multiqc -y
```

### 3. Clone this repository
```bash
git clone https://github.com/rarukua/sRNA-piRNA-pipeline.git
cd sRNA-piRNA-pipeline
```

---

## Reference Index Setup

### STAR genome index
```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles /path/to/hg38.fa \
     --sjdbGTFfile /path/to/annotation.gtf \
     --runThreadN 8
```

### Bowtie1 index
```bash
bowtie-build --threads 8 /path/to/hg38.fa /path/to/bowtie1_index/hg38
# Produces: hg38.1.ebwt  hg38.2.ebwt  hg38.3.ebwt  hg38.4.ebwt  hg38.rev.1.ebwt  hg38.rev.2.ebwt
```

### Bowtie2 index
```bash
bowtie2-build --threads 8 /path/to/hg38.fa /path/to/bowtie2_index/hg38
# Produces: hg38.1.bt2  hg38.2.bt2  hg38.3.bt2  hg38.4.bt2  hg38.rev.1.bt2  hg38.rev.2.bt2
```

---

## Usage

### Prepare samplesheet

```csv
sample_id,fastq
sample1,/absolute/path/to/sample1.fastq.gz
sample2,/absolute/path/to/sample2.fastq.gz
```

### Run with STAR (recommended)
```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/star_index \
    --gtf /path/to/piRbase_cluster.gtf \
    --aligner star \
    --outdir results
```

### Run with Bowtie1 (piRNA community standard)
```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --bowtie1_index /path/to/bowtie1_index/hg38 \
    --gtf /path/to/piRbase_cluster.gtf \
    --aligner bowtie1 \
    --outdir results
```

### Run with Bowtie2
```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --bowtie2_index /path/to/bowtie2_index/hg38 \
    --gtf /path/to/piRbase_cluster.gtf \
    --aligner bowtie2 \
    --outdir results
```

### Run with Docker
```bash
nextflow run main.nf -profile docker \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/star_index \
    --gtf /path/to/piRbase_cluster.gtf \
    --outdir results
```

### Run with Singularity (HPC)
```bash
nextflow run main.nf -profile singularity \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/star_index \
    --gtf /path/to/piRbase_cluster.gtf \
    --outdir results
```

---

## Output Structure

```
results/
├── fastqc/                        # Raw read QC
├── fastp/
│   ├── *.fastp.fastq.gz           # Trimmed reads
│   └── reports/                   # fastp HTML + JSON per sample
├── star/                          # BAM files + STAR logs (--aligner star)
├── bowtie1/                       # BAM files + logs (--aligner bowtie1)
├── bowtie2/                       # BAM files + logs (--aligner bowtie2)
├── featurecounts/                 # Count tables per sample
├── multiqc/
│   └── *_multiqc_report.html      # Aggregated QC report — open in browser
├── pipeline_report.html           # Nextflow execution report
├── pipeline_timeline.html         # Per-process timing
└── pipeline_dag.html              # Interactive workflow diagram
```

---

## Customising Parameters

All parameters can be overridden at runtime:

```bash
# Change piRNA size cutoff
nextflow run main.nf --max_len 31 ...

# Change minimum overlap for featureCounts
nextflow run main.nf --min_overlap 18 ...

# Change featureCounts feature type
nextflow run main.nf --fc_feature_type gene ...

# Run in QC-only mode (stop after fastp, before alignment)
nextflow run main.nf --mode qc_only ...
```

---

## piRNA Annotation

This pipeline uses **piRBase cluster GTF** format. Download from:
- [piRBase](http://www.pirbase.org/) — comprehensive piRNA database
- [piRNAdb](https://www.pirnai.org/) — alternative resource

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17).
- **STAR**: Dobin et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1).
- **Bowtie**: Langmead et al. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, 10(3).
- **Bowtie2**: Langmead & Salzberg (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4).
- **featureCounts**: Liao et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7).
- **MultiQC**: Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools and samples. *Bioinformatics*, 32(19).
- **Nextflow**: Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4).

---

## Author

**zgao** — PhD, Bioinformatics

Pipeline developed for reproducible piRNA-seq data analysis (hg38, piRBase annotation).
