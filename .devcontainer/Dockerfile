FROM rocker/r-ver:4.5.3
LABEL org.opencontainers.image.source=\
    https://github.com/mirramockakid/Mechanisms-of-CDK4-6-resistance

# Useful baseline tooling for analysis workflows
# in containers.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
        curl \
        ca-certificates \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager'); \
    BiocManager::install(c( \
        'limma', 'dplyr', 'tidyr', \
        'ggplot2', 'GEOquery', 'Rcpp', \
        'edgeR', 'assertthat', 'DESeq2', \
        'fgsea', 'msigdbr', \
        'Biobase', 'SummarizedExperiment', \
        'ComplexHeatmap', 'circlize', \
        'ggrepel', 'patchwork', \
        'purrr', 'readr', 'stringr' \
    ))"

# create dev user
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} \
        -m -s /bin/bash ${USERNAME}


