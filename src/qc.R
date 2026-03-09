# get raw counts, feature annotations and metadata for GSE222367
id <- "GSE222367"
gse <- GEOquery::getGEO(id, GSEMatrix = TRUE)

# download data files for GSE222367
GEOquery::getGEOSuppFiles(id, makeDirectory = FALSE, baseDir = "data/")

# list the contents of the downloaded tar file
members <- utils::untar("data/GSE222367_RAW.tar", list = TRUE)

# extract names GSM* of the files in the tar archive
utils::untar("data/GSE222367_RAW.tar", exdir = "data/")
files <- list.files("data/", pattern = "htseqcount.txt.gz", full.names = TRUE)
names(files) <- stringr::str_extract(members, "GSM\\d+")

# returns tibble of rows (features) and columns (samples) with raw counts
# column 1 is feature ID as Symbol
counts <- lapply(files, function(x) {
  print(x)
  readr::read_tsv(x,
    col_names = c("gene_id", "count")
  )
}) |>
  setNames(names(files)) |>
  purrr::imap(~ dplyr::rename(.x, !!.y := count)) |>
  purrr::reduce(dplyr::full_join, by = "gene_id")

# format metadata
meta <- Biobase::pData(gse[[1]]) |>
  dplyr::select(
    "Sample" = geo_accession, "Cell line" = `cell line:ch1`,
    "Treatment" = `treatment:ch1`, title
  ) |>
  dplyr::mutate("Group" = paste(Treatment, `Cell line`, sep = " ")) |>
  dplyr::mutate("Type" = stringr::str_extract(Treatment, "\\w+-resistant$")) |>
  dplyr::mutate(Type = ifelse(is.na(Type), "Sensitive", Type)) |>
  dplyr::mutate(`Dose uM` = stringr::str_extract(Treatment, "\\d+\\.\\d+") |> as.numeric()) |>
  dplyr::mutate(Resistance = stringr::str_extract(Treatment, "(Palbociclib|Abemaciclib|Ribociclib)")) |>
  dplyr::mutate(Pool = stringr::str_extract(title, "\\d+$")) |>
  dplyr::select(-title)


# ensure group variable is correctly formatted as factor with levels in the desired order
meta$Group <- factor(meta$Group)
levels(meta$Group) <- make.names(levels(meta$Group))

# check consistency of sample names between counts and metadata
assertthat::are_equal(all(colnames(counts)[-1] == meta$Sample), TRUE) |>
  assertthat::assert_that(msg = "Sample names in counts and metadata do not match")

# initialize DGEList object from edgeR package with raw counts
dge <- edgeR::DGEList(counts = counts)

# filter out genes with low expression
keep <- edgeR::filterByExpr(dge, group = meta$Group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# adds normalization factors to the DGEList object using the TMM method
dge <- edgeR::calcNormFactors(dge, method = "TMM")

# rlog transformation of counts for PCA and visualization
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = dge$counts,
  colData   = meta,
  design    = ~1
)

dds <- DESeq2::estimateSizeFactors(dds)

rld <- DESeq2::rlog(dds, blind = FALSE)
mat <- SummarizedExperiment::assay(rld)

# PCA plots
plot_pca(mat, meta, color_var = "Cell line", width = 6, height = 5)
plot_pca(mat, meta, color_var = "Treatment", width = 6, height = 5)
plot_pca(mat, meta, color_var = "Group", width = 6, height = 5)
plot_pca(mat, meta, color_var = "Dose uM", width = 6, height = 5)
plot_pca(mat, meta, color_var = "Resistance", width = 6, height = 5)
plot_pca(mat, meta, color_var = "Pool", width = 6, height = 5)

# density plots
plot_density_counts(mat, meta, color_var = "Group", width = 6, height = 5, is_raw_counts = FALSE)
plot_density_counts(mat, meta, color_var = "Cell line", width = 6, height = 5, is_raw_counts = FALSE)
plot_density_counts(mat, meta, color_var = "Treatment", width = 6, height = 5, is_raw_counts = FALSE)
plot_density_counts(mat, meta, color_var = "Dose", width = 6, height = 5, is_raw_counts = FALSE)
plot_density_counts(mat, meta, color_var = "Resistance", width = 6, height = 5, is_raw_counts = FALSE)

plot_ma(mat)
plot_ma(dge$counts)
