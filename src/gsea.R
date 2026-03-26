library(fgsea)
library(msigdbr)

# Hallmark gene sets as a named list for fgsea
hallmarks <- msigdbr(species = "Homo sapiens", collection = "H")
gene_sets <- split(hallmarks$gene_symbol, hallmarks$gs_name)

dir.create("tables/gsea", showWarnings = FALSE, recursive = TRUE)

gsea_list <- lapply(names(results_list), function(contrast_name) {

  res <- results_list[[contrast_name]]

  # Rank metric: -log10(P) * sign(logFC)
  scores <- -log10(res$P.Value) * sign(res$logFC)
  names(scores) <- res$Gene

  # Remove NAs and sort descending
  scores <- sort(scores[!is.na(scores)], decreasing = TRUE)

  set.seed(42)
  gsea_res <- fgseaMultilevel(
    pathways   = gene_sets,
    stats      = scores,
    minSize    = 15,
    maxSize    = 500
  )

  gsea_res <- gsea_res[order(gsea_res$pval), ]

  # GSEA plots – one PNG per pathway, ordered by p-value, in a contrast subdir
  plot_gsea(
    stats         = scores,
    pathways      = gene_sets,
    pathway_names = gsea_res$pathway,
    out_dir       = file.path("figures/gsea", contrast_name)
  )

  # Flatten the leading-edge gene list to a semicolon-separated string for CSV export
  gsea_res$leadingEdge <- vapply(gsea_res$leadingEdge, paste, character(1), collapse = ";")

  write.csv(
    gsea_res,
    file      = paste0("tables/gsea/", contrast_name, "_gsea_hallmarks.csv"),
    row.names = FALSE
  )

  gsea_res
}) |> setNames(names(results_list))
