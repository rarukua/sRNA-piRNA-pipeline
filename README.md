# sRNA / piRNA Analysis Pipeline

A reproducible Nextflow pipeline for small RNA-seq analysis with a focus on **piRNA biology**. Supports both STAR and Bowtie2 alignment, fastp-based size selection, and piRNA cluster quantification with featureCounts.

---

## Pipeline Overview

```
Raw FASTQ
    │
    ▼
 FastQC          ← QC on raw reads
    │
    ▼
 Fastp           ← Adapter trimming + size selection (≤24 nt for piRNA)
    │
    ▼
 STAR            ← Splice-aware, multi-mapper aware alignment
   or                (0 mismatches, up to 1M multimappers, local mode)
 Bowtie2         ← Alternative short-read aligner
    │
    ▼
 featureCounts   ← piRNA cluster quantification
                    (-M -O --minOverlap 22 against piRBase GTF)
    │
    ▼
 MultiQC         ← Aggregated QC report
```

---

## Key piRNA-Specific Parameters

| Step | Parameter | Value | Rationale |
|------|-----------|-------|-----------|
| fastp | `--max_len1` | 24 | Retain reads ≤24 nt (piRNA size range) |
| STAR | `--outFilterMismatchNmax` | 0 | Exact matches only |
| STAR | `--outFilterMultimapNmax` | 1,000,000 | Allow multi-mappers (piRNAs map to repeats) |
| STAR | `--outFilterMatchNmin` | 22 | Minimum matched bases |
| STAR | `--alignEndsType` | Local | Soft-clipping allowed |
| featureCounts | `--minOverlap` | 22 | Minimum read-feature overlap |
| featureCounts | `-M -O` | enabled | Count multi-mappers + overlapping features |
| featureCounts | `-Q 0` | 0 | No MAPQ filter (needed for multimappers) |

---

## Requirements

### Option A: Docker (recommended for reproducibility)
- [Docker](https://docs.docker.com/get-docker/) ≥ 20.x
- [Nextflow](https://www.nextflow.io/) ≥ 23.x

### Option B: Singularity (recommended on shared HPC/Linux servers)
- [Singularity](https://sylabs.io/singularity/) ≥ 3.x
- [Nextflow](https://www.nextflow.io/) ≥ 23.x

### Option C: Local tools (all must be in PATH)
- FastQC ≥ 0.12
- fastp ≥ 0.23
- STAR ≥ 2.7.11
- Bowtie2 ≥ 2.5
- SAMtools ≥ 1.19
- featureCounts (Subread ≥ 2.0)
- MultiQC ≥ 1.21
- Nextflow ≥ 23.x

---

## Installation

### 1. Install Nextflow

```bash
# Requires Java 11+
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### 2. Clone this repository

```bash
git clone https://github.com/YOUR_USERNAME/sRNA-piRNA-pipeline.git
cd sRNA-piRNA-pipeline
```

### 3. Build Docker image (if using Docker profile)

```bash
docker build -t zgao/pirna-pipeline:latest .
```

---

## Usage

### Prepare your samplesheet

Create a CSV file (see `assets/samplesheet.csv` for template):

```csv
sample_id,fastq
sample1,/absolute/path/to/sample1.fastq.gz
sample2,/absolute/path/to/sample2.fastq.gz
```

> **Note:** Single-end reads only. piRNA-seq is almost always single-end.

### Run with STAR (default)

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/STAR_index \
    --gtf /path/to/piRbase_cluster_GTF.gtf \
    --aligner star \
    --outdir results
```

### Run with Bowtie2

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --bowtie2_index /path/to/bowtie2_index/genome \
    --gtf /path/to/piRbase_cluster_GTF.gtf \
    --aligner bowtie2 \
    --outdir results
```

### Run with Docker

```bash
nextflow run main.nf \
    -profile docker \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/STAR_index \
    --gtf /path/to/piRbase_cluster_GTF.gtf \
    --outdir results
```

### Run with Singularity (HPC/Linux server)

```bash
nextflow run main.nf \
    -profile singularity \
    --samplesheet assets/samplesheet.csv \
    --genome_dir /path/to/STAR_index \
    --gtf /path/to/piRbase_cluster_GTF.gtf \
    --outdir results
```

---

## Output Structure

```
results/
├── fastqc/                  # Raw read QC (HTML + zip per sample)
├── fastp/
│   ├── *.fastp.fastq.gz     # Trimmed reads
│   └── reports/             # fastp HTML + JSON reports
├── star/                    # BAM files + STAR logs (if --aligner star)
├── bowtie2/                 # BAM files + alignment logs (if --aligner bowtie2)
├── featurecounts/           # Count tables per sample
├── multiqc/
│   ├── multiqc_report.html  # Aggregated QC report (open in browser)
│   └── multiqc_data/
├── pipeline_report.html     # Nextflow execution report
├── pipeline_timeline.html   # Per-process timing
└── pipeline_dag.svg         # Workflow DAG visualization
```

---

## Customising Parameters

All parameters can be overridden at runtime with `--param value`:

```bash
# Change piRNA size cutoff (e.g. include 26-31 nt piRNAs)
nextflow run main.nf --max_len 31 ...

# Change minimum overlap for featureCounts
nextflow run main.nf --min_overlap 18 ...

# Change featureCounts feature type (if your GTF uses different 3rd column)
nextflow run main.nf --fc_feature_type gene ...
```

---

## Reference Files

### STAR genome index

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/STAR_index \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/annotation.gtf \
     --genomeSAindexNbases 11    # use 11 for small genomes
```

### piRNA annotation

This pipeline uses **piRBase cluster GTF** format. Download from:
- [piRBase](http://www.pirbase.org/) — comprehensive piRNA database
- [piRNAdb](https://www.pirnai.org/) — alternative resource

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17).
- **STAR**: Dobin et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1).
- **Bowtie2**: Langmead & Salzberg (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4).
- **featureCounts**: Liao et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7).
- **MultiQC**: Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools and samples. *Bioinformatics*, 32(19).
- **Nextflow**: Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4).

---

## Author

**zgao** — PhD, Bioinformatics

Pipeline developed for reproducible piRNA-seq data analysis.
