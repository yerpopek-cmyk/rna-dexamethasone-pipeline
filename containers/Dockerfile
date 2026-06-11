# =============================================================================
# Dockerfile — rna-dexamethasone pipeline
# Base: rocker/r-ver:4.5 + bioinformatics tools
# =============================================================================

FROM rocker/r-ver:4.5

LABEL maintainer="your@email.com"
LABEL description="RNA-seq pipeline for dexamethasone effects analysis (GSE52778)"
LABEL version="1.0.0"

# ── System dependencies ───────────────────────────────────────────────────────
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    bzip2 \
    gzip \
    pigz \
    samtools \
    openjdk-17-jre-headless \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── FastQC ───────────────────────────────────────────────────────────────────
ARG FASTQC_VERSION=0.12.1
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
    -O /tmp/fastqc.zip \
    && unzip /tmp/fastqc.zip -d /opt/ \
    && chmod +x /opt/FastQC/fastqc \
    && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
    && rm /tmp/fastqc.zip

# ── fastp ────────────────────────────────────────────────────────────────────
ARG FASTP_VERSION=0.23.4
RUN wget -q https://github.com/OpenGene/fastp/releases/download/v${FASTP_VERSION}/fastp \
    -O /usr/local/bin/fastp \
    && chmod +x /usr/local/bin/fastp

# ── STAR aligner ─────────────────────────────────────────────────────────────
ARG STAR_VERSION=2.7.11b
RUN wget -q https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip \
    -O /tmp/STAR.zip \
    && unzip /tmp/STAR.zip -d /tmp/ \
    && cp /tmp/STAR_${STAR_VERSION}/Linux_x86_64_static/STAR /usr/local/bin/ \
    && chmod +x /usr/local/bin/STAR \
    && rm -rf /tmp/STAR*

# ── Subread (featureCounts) ───────────────────────────────────────────────────
ARG SUBREAD_VERSION=2.0.6
RUN wget -q https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-source.tar.gz \
    -O /tmp/subread-source.tar.gz \
    && tar xzf /tmp/subread-source.tar.gz -C /opt/ \
    && cd /opt/subread-${SUBREAD_VERSION}-source/src \
    && make -f Makefile.Linux \
    && ln -s /opt/subread-${SUBREAD_VERSION}-source/bin/featureCounts /usr/local/bin/featureCounts \
    && rm /tmp/subread-source.tar.gz

# ── MultiQC ──────────────────────────────────────────────────────────────────
RUN pip3 install --no-cache-dir multiqc==1.21

# ── R CRAN packages ──────────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
    'optparse', \
    'ggplot2', \
    'dplyr', \
    'tibble', \
    'tidyr', \
    'readr', \
    'scales', \
    'RColorBrewer', \
    'pheatmap', \
    'rmarkdown', \
    'knitr', \
    'kableExtra', \
    'DT', \
    'cowplot', \
    'ggrepel', \
    'stringr', \
    'viridis' \
    ), repos='https://cloud.r-project.org', Ncpus=4)"

# ── R Bioconductor packages ───────────────────────────────────────────────────
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'ReactomePA', \
    'clusterProfiler', \
    'enrichplot', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'EnhancedVolcano', \
    'AnnotationDbi', \
    'biomaRt', \
    'fgsea', \
    'DOSE', \
    'ggupset' \
    ), update=FALSE, ask=FALSE)"

# Install perl and pandoc at the end to keep R packages layer cache intact
RUN apt-get update && apt-get install -y --no-install-recommends perl pandoc && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Set working directory ─────────────────────────────────────────────────────
WORKDIR /data

# ── Default command ───────────────────────────────────────────────────────────
CMD ["bash"]
