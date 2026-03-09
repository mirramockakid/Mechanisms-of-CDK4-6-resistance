# set library path to a local directory
dir.create("rlib", showWarnings = FALSE)
.libPaths(normalizePath("./rlib"))

# setup folder structure
dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "rlib")
}

BiocManager::install("limma", lib = "rlib")
BiocManager::install("tidyverse", lib = "rlib")
BiocManager::install("ggplot2", lib = "rlib")
BiocManager::install("GEOquery", lib = "rlib")
BiocManager::install("Rcpp", lib = "rlib")
BiocManager::install("edgeR", lib = "rlib")
BiocManager::install("assertthat", lib = "rlib")
BiocManager::install("DESeq2", lib = "rlib")

source("src/plotting.R")
