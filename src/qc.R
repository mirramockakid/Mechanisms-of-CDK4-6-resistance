# get raw counts, feature annotations and metadata for GSE222367
id <- "GSE222367"
gse <- GEOquery::getGEO(id, GSEMatrix = TRUE)
raw_counts <- exprs(gse[[1]])
feature_annotations <- fData(gse[[1]])
metadata <- pData(gse[[1]])

# download data files for GSE222367
GEOquery::getGEOSuppFiles(id, makeDirectory = FALSE, baseDir = "data/")

# list the contents of the downloaded tar file
members <- utils::untar("data/GSE222367_RAW.tar", list = TRUE)

# extract names GSM* of the files in the tar archive
utils::untar("data/GSE222367_RAW.tar", exdir = "data/"))
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
  setNames(members) |>
  purrr::imap(~ dplyr::rename(.x, !!.y := count)) |>
  purrr::reduce(dplyr::full_join, by = "gene_id")

# format metadata
x <- pData(gse[[1]]) |>
  dplyr::select("Sample" = geo_accession, "Cell line" = `cell line:ch1`, "Group" = `treatment:ch1`)


x$Group
