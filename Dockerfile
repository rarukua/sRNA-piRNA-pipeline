# ============================================================
# Dockerfile — sRNA / piRNA Analysis Pipeline
# Tools: FastQC, fastp, STAR, Bowtie2, SAMtools,
#        Subread (featureCounts), MultiQC
# Base:  Ubuntu 22.04 (LTS, widely supported on HPC)
# ============================================================

FROM ubuntu:22.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/conda/bin:${PATH}"

LABEL maintainer="zgao"
LABEL description="piRNA small RNA-seq analysis pipeline"
LABEL version="1.0"

# --- System dependencies ---
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    unzip \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# --- Install Miniconda (makes tool versioning reproducible) ---
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && conda clean -afy

# --- Configure conda channels ---
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# --- Install all bioinformatics tools in one layer ---
# Pinning versions for reproducibility
RUN conda install -y \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    star=2.7.11a \
    bowtie2=2.5.3 \
    samtools=1.19 \
    subread=2.0.6 \
    multiqc=1.21 \
    && conda clean -afy

# --- Set working directory ---
WORKDIR /data

# --- Default command: show help ---
CMD ["bash"]

# ============================================================
# Build instructions:
#   docker build -t zgao/pirna-pipeline:latest .
#
# Run example:
#   docker run --rm -v $(pwd):/data zgao/pirna-pipeline:latest \
#       nextflow run /data/main.nf -profile docker
# ============================================================
