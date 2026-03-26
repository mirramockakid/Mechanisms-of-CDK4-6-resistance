plot_pca <- function(expr_mat, metadata, color_var,
                     out_dir = "figures/pca",
                     width = 7, height = 5,
                     normalised = FALSE) {
  if (!is.matrix(expr_mat)) {
    stop("`expr_mat` must be a matrix (e.g., rlog assay matrix).")
  }

  if (!is.numeric(expr_mat)) {
    stop("`expr_mat` must be a numeric matrix.")
  }

  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data.frame.")
  }

  if (!is.character(color_var) || length(color_var) != 1) {
    stop("`color_var` must be a single character string.")
  }

  if (!is.numeric(width) || length(width) != 1 || width <= 0) {
    stop("`width` must be a single positive number.")
  }

  if (!is.numeric(height) || length(height) != 1 || height <= 0) {
    stop("`height` must be a single positive number.")
  }

  if (!color_var %in% colnames(metadata)) {
    stop(sprintf("`%s` was not found in metadata columns.", color_var))
  }

  sample_ids <- colnames(expr_mat)
  if (is.null(sample_ids)) {
    stop("`expr_mat` must have column names corresponding to sample IDs.")
  }

  if (!"Sample" %in% colnames(metadata)) {
    stop("`metadata` must contain a `Sample` column",
         " matching expression matrix sample IDs.")
  }

  if (!all(sample_ids %in% metadata$Sample)) {
    stop("Not all expression matrix sample IDs were found in metadata$Sample.")
  }

  metadata <- metadata[match(sample_ids, metadata$Sample), , drop = FALSE]
  pca <- stats::prcomp(t(expr_mat), center = TRUE, scale. = TRUE)

  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pca_df <- dplyr::left_join(pca_df, metadata, by = "Sample")
  color_data <- pca_df[[color_var]]

  pc_var <- (pca$sdev^2) / sum(pca$sdev^2)
  xlab <- sprintf("PC1 (%.1f%%)", 100 * pc_var[1])
  ylab <- sprintf("PC2 (%.1f%%)", 100 * pc_var[2])

  plot_obj <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = .data[[color_var]]
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "grey85") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey85") +
    ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::labs(
      title = "Principal Component Analysis",
      subtitle = sprintf("Points colored by %s", color_var),
      x = xlab,
      y = ylab,
      color = color_var
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (is.numeric(color_data)) {
    plot_obj <- plot_obj +
      ggplot2::scale_color_gradientn(
        colors = grDevices::hcl.colors(9, "TealGrn")
      )
  } else {
    n_levels <- length(unique(as.character(color_data)))
    plot_obj <- plot_obj +
      ggplot2::scale_color_manual(
        values = grDevices::hcl.colors(n_levels, "Dark 3")
      )
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", color_var)
  if (normalised) {
    safe_name <- paste0(safe_name, "_norm")
  }
  out_file <- file.path(out_dir, sprintf("%s.pdf", safe_name))
  ggplot2::ggsave(
    filename = out_file,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height
  )

  invisible(plot_obj)
}

plot_density_counts <- function(expr_mat, metadata,
                               color_var,
                               out_dir = "figures/density",
                               width = 8, height = 6,
                               is_raw_counts = NULL,
                               normalised = FALSE) {
  if (!is.matrix(expr_mat)) {
    stop("`expr_mat` must be a matrix.")
  }

  if (!is.numeric(expr_mat)) {
    stop("`expr_mat` must be a numeric matrix.")
  }

  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data.frame.")
  }

  if (!is.character(color_var) || length(color_var) != 1) {
    stop("`color_var` must be a single character string.")
  }

  if (!is.numeric(width) || length(width) != 1 || width <= 0) {
    stop("`width` must be a single positive number.")
  }

  if (!is.numeric(height) || length(height) != 1 || height <= 0) {
    stop("`height` must be a single positive number.")
  }

  if (!is.null(is_raw_counts) &&
      (!is.logical(is_raw_counts) ||
       length(is_raw_counts) != 1)) {
    stop("`is_raw_counts` must be NULL or a single TRUE/FALSE value.")
  }

  if (!color_var %in% colnames(metadata)) {
    stop(sprintf("`%s` was not found in metadata columns.", color_var))
  }

  sample_ids <- colnames(expr_mat)
  if (is.null(sample_ids)) {
    stop("`expr_mat` must have column names corresponding to sample IDs.")
  }

  if (!"Sample" %in% colnames(metadata)) {
    stop("`metadata` must contain a `Sample` column",
         " matching count matrix sample IDs.")
  }

  if (!all(sample_ids %in% metadata$Sample)) {
    stop("Not all count matrix sample IDs were found in metadata$Sample.")
  }

  metadata <- metadata[match(sample_ids, metadata$Sample), , drop = FALSE]
  if (is.null(is_raw_counts)) {
    finite_vals <- expr_mat[is.finite(expr_mat)]
    is_raw_counts <- length(finite_vals) > 0 &&
      all(finite_vals >= 0) &&
      all(abs(finite_vals - round(finite_vals)) < .Machine$double.eps^0.5)
  }

  plot_mat <- if (is_raw_counts) log2(expr_mat + 1) else expr_mat
  x_label <- if (is_raw_counts) "log2(count + 1)" else "Expression"

  density_df <- do.call(
    rbind,
    lapply(seq_along(sample_ids), function(i) {
      data.frame(
        Sample = sample_ids[i],
        Expression = plot_mat[, i],
        stringsAsFactors = FALSE
      )
    })
  )

  density_df <- dplyr::left_join(density_df, metadata, by = "Sample")
  color_data <- density_df[[color_var]]

  plot_obj <- ggplot2::ggplot(
    density_df,
    ggplot2::aes(
      x = Expression,
      group = Sample,
      color = .data[[color_var]]
    )
  ) +
    ggplot2::geom_density(linewidth = 0.6, alpha = 0.9) +
    ggplot2::labs(
      title = "Sample Expression Density",
      subtitle = sprintf("Density curves colored by %s", color_var),
      x = x_label,
      y = "Density",
      color = color_var
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (is.numeric(color_data)) {
    plot_obj <- plot_obj +
      ggplot2::scale_color_gradientn(
        colors = grDevices::hcl.colors(9, "TealGrn")
      )
  } else {
    n_levels <- length(unique(as.character(color_data)))
    plot_obj <- plot_obj +
      ggplot2::scale_color_manual(
        values = grDevices::hcl.colors(n_levels, "Dark 3")
      )
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", color_var)
  if (normalised) {
    safe_name <- paste0(safe_name, "_norm")
  }
  out_file <- file.path(out_dir, sprintf("%s.pdf", safe_name))
  ggplot2::ggsave(
    filename = out_file,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height
  )

  invisible(plot_obj)
}

plot_ma <- function(counts_mat, norm_mat = NULL,
                   out_dir = "figures/ma",
                   file_name = "ma_plots.pdf",
                   plots_per_page = 10, ncol = 2,
                   width = 11, height = 14) {
  if (!is.matrix(counts_mat)) {
    stop("`counts_mat` must be a matrix.")
  }

  if (!is.numeric(counts_mat)) {
    stop("`counts_mat` must be a numeric matrix.")
  }

  if (!is.null(norm_mat)) {
    if (!is.matrix(norm_mat)) {
      stop("`norm_mat` must be a matrix.")
    }
    if (!is.numeric(norm_mat)) {
      stop("`norm_mat` must be a numeric matrix.")
    }
    if (!identical(dim(counts_mat), dim(norm_mat))) {
      stop("`counts_mat` and `norm_mat` must have the same dimensions.")
    }
    if (!identical(colnames(counts_mat), colnames(norm_mat))) {
      stop("`counts_mat` and `norm_mat` must have the same column names.")
    }
  }

  if (!is.numeric(width) || length(width) != 1 || width <= 0) {
    stop("`width` must be a single positive number.")
  }

  if (!is.numeric(height) || length(height) != 1 || height <= 0) {
    stop("`height` must be a single positive number.")
  }

  if (!is.character(file_name) || length(file_name) != 1) {
    stop("`file_name` must be a single character string.")
  }

  if (!is.numeric(plots_per_page) ||
      length(plots_per_page) != 1 ||
      plots_per_page <= 0) {
    stop("`plots_per_page` must be a single positive number.")
  }

  if (!is.numeric(ncol) || length(ncol) != 1 || ncol <= 0) {
    stop("`ncol` must be a single positive number.")
  }

  sample_ids <- colnames(counts_mat)
  if (is.null(sample_ids)) {
    stop("`counts_mat` must have column names corresponding to sample IDs.")
  }

  log_counts <- log2(counts_mat + 1)
  ref_counts <- rowMeans(log_counts, na.rm = TRUE)

  if (!is.null(norm_mat)) {
    log_norm <- norm_mat
    ref_norm <- rowMeans(log_norm, na.rm = TRUE)
  } else {
    log_norm <- log_counts
    ref_norm <- ref_counts
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  m_limit <- 5

  plot_list <- vector("list", length(sample_ids))
  names(plot_list) <- sample_ids

  for (i in seq_along(sample_ids)) {
    sample_name <- sample_ids[i]

    plot_df <- data.frame(
      A = 0.5 * (log_counts[, i] + ref_counts),
      M = log_norm[, i] - ref_norm,
      stringsAsFactors = FALSE
    )
    plot_df <- plot_df[
      is.finite(plot_df$A) & is.finite(plot_df$M), ,
      drop = FALSE
    ]
    plot_df$M <- pmax(pmin(plot_df$M, m_limit), -m_limit)

    plot_obj <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = A, y = M)
    ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "grey85") +
      ggplot2::geom_point(size = 1.1, alpha = 0.55, color = "#0072B2") +
      ggplot2::geom_smooth(
        method = "loess", color = "red",
        linewidth = 0.6, se = FALSE
      ) +
      ggplot2::labs(
        title = sample_name,
        x = "A",
        y = "M"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          color = "grey90", linewidth = 0.3
        ),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::coord_cartesian(ylim = c(-m_limit, m_limit))

    plot_list[[i]] <- plot_obj
  }


  plots_per_page <- as.integer(plots_per_page)
  ncol <- as.integer(ncol)
  nrow <- ceiling(plots_per_page / ncol)
  safe_file <- gsub("[^A-Za-z0-9._-]", "_", file_name)
  if (!grepl("\\.pdf$", safe_file, ignore.case = TRUE)) {
    safe_file <- paste0(safe_file, ".pdf")
  }
  out_file <- file.path(out_dir, safe_file)

  grDevices::pdf(out_file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (start_idx in seq(1, length(plot_list), by = plots_per_page)) {
    end_idx <- min(start_idx + plots_per_page - 1, length(plot_list))
    page_plots <- plot_list[start_idx:end_idx]

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(nrow = nrow, ncol = ncol)
    ))

    for (j in seq_along(page_plots)) {
      row_idx <- ((j - 1) %/% ncol) + 1
      col_idx <- ((j - 1) %% ncol) + 1
      print(
        page_plots[[j]],
        vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = col_idx)
      )
    }
  }

  invisible(plot_list)
}

plot_volcano <- function(df,
                         out_dir = "figures/volcano",
                         file_name = "volcano.pdf",
                         title = "Volcano Plot",
                         width = 7, height = 6) {
  required <- c("Gene", "logFC", "P.Value", "adj.P.Val")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(
      "Missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }

  if (!is.numeric(width) ||
      length(width) != 1 || width <= 0) {
    stop("`width` must be a single positive number.")
  }

  if (!is.numeric(height) ||
      length(height) != 1 || height <= 0) {
    stop("`height` must be a single positive number.")
  }

  if (!is.character(file_name) ||
      length(file_name) != 1) {
    stop("`file_name` must be a single character string.")
  }

  df$neg_log10_p <- -log10(df$P.Value)

  adj_cutoff <- max(
    df$P.Value[df$adj.P.Val < 0.05],
    na.rm = TRUE
  )
  hline_y <- -log10(adj_cutoff)

  plot_obj <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = logFC, y = neg_log10_p)
  ) +
    ggplot2::geom_point(
      size = 1.2, alpha = 0.5, color = "grey40"
    ) +
    ggplot2::geom_hline(
      yintercept = hline_y,
      linewidth = 0.5,
      linetype = "dashed",
      color = "red"
    ) +
    ggrepel::geom_text_repel(
      data = utils::head(
        df[order(-df$neg_log10_p), ], 20
      ),
      ggplot2::aes(label = Gene),
      size = 3, max.overlaps = 20
    ) +
    ggplot2::labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(P-value)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      )
    )

  dir.create(
    out_dir, showWarnings = FALSE, recursive = TRUE
  )
  safe_file <- gsub(
    "[^A-Za-z0-9._-]", "_", file_name
  )
  if (!grepl("\\.pdf$", safe_file, ignore.case = TRUE)) {
    safe_file <- paste0(safe_file, ".pdf")
  }
  out_file <- file.path(out_dir, safe_file)
  ggplot2::ggsave(
    filename = out_file,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height
  )

  invisible(plot_obj)
}

plot_heatmap <- function(expr_mat, metadata, results,
                         out_dir = "figures/heatmap",
                         file_name = "heatmap.pdf",
                         pval_threshold = 0.05,
                         max_genes = 50,
                         annotation = NULL,
                         gene_names = NULL,
                         width = 10, height = 8) {
  if (!is.matrix(expr_mat) || !is.numeric(expr_mat))
    stop("`expr_mat` must be a numeric matrix.")
  if (!is.data.frame(metadata))
    stop("`metadata` must be a data.frame.")
  if (!is.data.frame(results))
    stop("`results` must be a data.frame.")

  required_cols <- c("Gene", "adj.P.Val", "logFC")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0)
    stop("Missing columns in `results`: ", paste(missing_cols, collapse = ", "))

  if (!"Sample" %in% colnames(metadata))
    stop("`metadata` must contain a `Sample` column.")

  if (!is.null(annotation)) {
    bad <- setdiff(annotation, colnames(metadata))
    if (length(bad) > 0)
      stop("Annotation column(s) not found in metadata: ", paste(bad, collapse = ", "))
  }

  if (!is.null(gene_names)) {
    if (length(gene_names) != nrow(expr_mat))
      stop("`gene_names` must have the same length as the number of rows in `expr_mat`.")
    rownames(expr_mat) <- gene_names
  }

  # Select significant genes
  sig <- results[!is.na(results$adj.P.Val) & results$adj.P.Val < pval_threshold, ]
  if (nrow(sig) == 0) {
    message("No significant genes at adj.P.Val < ", pval_threshold, "; skipping heatmap.")
    return(invisible(NULL))
  }
  if (nrow(sig) > max_genes)
    sig <- sig[order(sig$adj.P.Val), ][seq_len(max_genes), ]

  genes <- sig$Gene
  mat <- expr_mat[rownames(expr_mat) %in% genes, , drop = FALSE]

  # Row-wise z-score
  mat_z <- t(scale(t(mat)))

  # Column annotation
  col_ann <- NULL
  if (!is.null(annotation)) {
    ann_df <- metadata[match(colnames(mat_z), metadata$Sample), annotation, drop = FALSE]
    rownames(ann_df) <- colnames(mat_z)

    col_list <- stats::setNames(lapply(annotation, function(v) {
      vals <- ann_df[[v]]
      if (is.numeric(vals)) {
        circlize::colorRamp2(range(vals, na.rm = TRUE), c("white", "#08306B"))
      } else {
        lvls <- unique(na.omit(as.character(vals)))
        stats::setNames(grDevices::hcl.colors(length(lvls), "Dark 3"), lvls)
      }
    }), annotation)

    ha_data <- stats::setNames(lapply(annotation, function(v) {
      vals <- ann_df[[v]]
      if (is.numeric(vals)) vals else as.character(vals)
    }), annotation)

    col_ann <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(ha_data, list(col = col_list, annotation_name_side = "left"))
    )
  }

  lim <- max(abs(mat_z), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))

  ht <- ComplexHeatmap::Heatmap(
    mat_z,
    name = "Z-score",
    col = col_fun,
    top_annotation = col_ann,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8)
  )

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_file <- gsub("[^A-Za-z0-9._-]", "_", file_name)
  if (!grepl("\\.pdf$", safe_file, ignore.case = TRUE))
    safe_file <- paste0(safe_file, ".pdf")
  out_file <- file.path(out_dir, safe_file)

  grDevices::pdf(out_file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)

  invisible(ht)
}

plot_scatter <- function(res1, res2,
                         name1, name2,
                         out_dir = "figures/scatter",
                         file_name = "scatter.pdf",
                         label_top_n = 10,
                         width = 7, height = 6) {
  required_cols <- c("Gene", "logFC", "adj.P.Val")

  missing1 <- setdiff(required_cols, colnames(res1))
  if (length(missing1) > 0)
    stop("Missing columns in `res1`: ", paste(missing1, collapse = ", "))

  missing2 <- setdiff(required_cols, colnames(res2))
  if (length(missing2) > 0)
    stop("Missing columns in `res2`: ", paste(missing2, collapse = ", "))

  if (!is.character(name1) || length(name1) != 1)
    stop("`name1` must be a single character string.")

  if (!is.character(name2) || length(name2) != 1)
    stop("`name2` must be a single character string.")

  if (!is.numeric(width) || length(width) != 1 || width <= 0)
    stop("`width` must be a single positive number.")

  if (!is.numeric(height) || length(height) != 1 || height <= 0)
    stop("`height` must be a single positive number.")

  if (!is.numeric(label_top_n) || length(label_top_n) != 1 || label_top_n < 0)
    stop("`label_top_n` must be a single non-negative number.")

  merged <- merge(
    res1[, c("Gene", "logFC", "adj.P.Val")],
    res2[, c("Gene", "logFC", "adj.P.Val")],
    by = "Gene",
    suffixes = c(".c1", ".c2")
  )

  sig1 <- !is.na(merged$adj.P.Val.c1) & merged$adj.P.Val.c1 < 0.05
  sig2 <- !is.na(merged$adj.P.Val.c2) & merged$adj.P.Val.c2 < 0.05

  merged$significance <- dplyr::case_when(
    sig1 & sig2  ~ "Both",
    sig1 & !sig2 ~ "Contrast 1 only",
    !sig1 & sig2 ~ "Contrast 2 only",
    TRUE         ~ "Not significant"
  )
  merged$significance <- factor(
    merged$significance,
    levels = c("Both", "Contrast 1 only", "Contrast 2 only", "Not significant")
  )

  # Alpha scaled to the mean absolute logFC across the two contrasts
  merged$mean_abs_logFC <- (abs(merged$logFC.c1) + abs(merged$logFC.c2)) / 2
  max_abs <- max(merged$mean_abs_logFC, na.rm = TRUE)
  merged$pt_alpha <- 0.1 + 0.9 * (merged$mean_abs_logFC / max_abs)

  # Label top/bottom n genes on each axis (up to 4 * label_top_n unique genes)
  n <- as.integer(label_top_n)
  label_genes <- unique(c(
    utils::head(merged$Gene[order(-merged$logFC.c1)], n),  # top x
    utils::head(merged$Gene[order( merged$logFC.c1)], n),  # bottom x
    utils::head(merged$Gene[order(-merged$logFC.c2)], n),  # top y
    utils::head(merged$Gene[order( merged$logFC.c2)], n)   # bottom y
  ))
  label_df <- merged[merged$Gene %in% label_genes, ]

  colour_map <- c(
    "Both"             = "red",
    "Contrast 1 only"  = "green",
    "Contrast 2 only"  = "blue",
    "Not significant"  = "grey70"
  )

  plot_obj <- ggplot2::ggplot(
    merged,
    ggplot2::aes(x = logFC.c1, y = logFC.c2, colour = significance)
  ) +
    ggplot2::geom_point(size = 1.2, alpha = merged$pt_alpha) +
    ggplot2::geom_smooth(
      method      = "lm",
      formula     = y ~ x,
      colour      = "black",
      linewidth   = 0.7,
      linetype    = "dotted",
      se          = FALSE,
      inherit.aes = FALSE,
      ggplot2::aes(x = logFC.c1, y = logFC.c2)
    ) +
    ggplot2::scale_colour_manual(
      values = colour_map,
      name = "Significance (adj.P.Val < 0.05)"
    ) +
    ggplot2::labs(
      x = sprintf("logFC (%s)", name1),
      y = sprintf("logFC (%s)", name2)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.title      = ggplot2::element_text(face = "bold"),
      axis.title        = ggplot2::element_text(face = "bold")
    )

  if (nrow(label_df) > 0) {
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(
        data          = label_df,
        ggplot2::aes(label = Gene),
        colour        = "black",
        size          = 3,
        max.overlaps  = 20,
        show.legend   = FALSE
      )
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_file <- gsub("[^A-Za-z0-9._-]", "_", file_name)
  if (!grepl("\\.pdf$", safe_file, ignore.case = TRUE))
    safe_file <- paste0(safe_file, ".pdf")
  out_file <- file.path(out_dir, safe_file)
  ggplot2::ggsave(
    filename = out_file,
    plot     = plot_obj,
    device   = "pdf",
    width    = width,
    height   = height
  )

  invisible(plot_obj)
}

# ---------------------------------------------------------------------------
# plot_gsea
#
# Produces a classic three-panel GSEA running-enrichment-score plot for one
# or more pathways.  Each pathway is rendered on its own page of the output
# PDF.
#
# Arguments:
#   stats         Named, sorted (decreasing) numeric rank-metric vector.
#   pathways      Named list of character vectors (gene sets), OR a single
#                 character vector for one pathway.
#   pathway_names Character vector of names to plot.  When NULL (default),
#                 all entries in `pathways` are plotted.
#   out_dir       Output directory (created if absent).
#   file_name     Output file name (must end in .pdf; extension added if
#                 missing).
#   width, height Page dimensions in inches.
# ---------------------------------------------------------------------------
plot_gsea <- function(stats,
                      pathways,
                      pathway_names = NULL,
                      out_dir       = "figures/gsea",
                      width         = 7,
                      height        = 7) {

  # ---- input validation ---------------------------------------------------
  if (!is.numeric(stats) || is.null(names(stats)))
    stop("`stats` must be a named numeric vector.")

  if (is.character(pathways))
    pathways <- list(pathway = pathways)

  if (!is.list(pathways) || length(pathways) == 0)
    stop("`pathways` must be a non-empty named list of character vectors.")

  if (is.null(names(pathways)))
    names(pathways) <- paste0("pathway_", seq_along(pathways))

  if (!is.null(pathway_names)) {
    missing_pw <- setdiff(pathway_names, names(pathways))
    if (length(missing_pw) > 0)
      stop("Pathway(s) not found in `pathways`: ",
           paste(missing_pw, collapse = ", "))
    pathways <- pathways[pathway_names]
  }

  if (!is.numeric(width)  || length(width)  != 1 || width  <= 0)
    stop("`width` must be a single positive number.")
  if (!is.numeric(height) || length(height) != 1 || height <= 0)
    stop("`height` must be a single positive number.")

  # ---- helpers -------------------------------------------------------------

  # Compute running enrichment score using the classical weighted KS statistic
  # (abs(stat)-weighted hits, uniform misses).
  .running_es <- function(stats, gene_set) {
    N     <- length(stats)
    in_set <- names(stats) %in% gene_set
    k      <- sum(in_set)

    if (k == 0L || k == N)
      return(rep(0, N))

    stat_sum <- sum(abs(stats[in_set]))
    if (stat_sum == 0) stat_sum <- 1   # guard against all-zero stats

    step_hit  <-  abs(stats) / stat_sum   # added when gene is in set
    step_miss <- -1 / (N - k)             # subtracted when gene is not in set

    increments <- ifelse(in_set, step_hit, step_miss)
    cumsum(increments)
  }

  # Build the three-panel ggplot for one pathway
  .one_pathway_plot <- function(pathway_name, gene_set, stats) {
    N      <- length(stats)
    in_set <- names(stats) %in% gene_set
    k      <- sum(in_set)

    running <- .running_es(stats, gene_set)
    es      <- running[which.max(abs(running))]
    es_pos  <- which.max(abs(running))

    # Enrichment score sign determines colour
    es_col <- if (es >= 0) "#B2182B" else "#2166AC"

    # --- Panel 1: running enrichment score ----------------------------------
    es_df <- data.frame(x = seq_len(N), es = running)

    p_es <- ggplot2::ggplot(es_df, ggplot2::aes(x = x, y = es)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey70") +
      ggplot2::geom_line(colour = es_col, linewidth = 0.9) +
      ggplot2::geom_vline(
        xintercept = es_pos,
        linewidth  = 0.5,
        linetype   = "dashed",
        colour     = "grey40"
      ) +
      ggplot2::annotate(
        "text",
        x      = es_pos,
        y      = es,
        label  = sprintf("ES = %.3f", es),
        hjust  = if (es_pos > N / 2) 1.1 else -0.1,
        vjust  = if (es >= 0) 1.5 else -0.5,
        size   = 3.5,
        colour = es_col
      ) +
      ggplot2::labs(
        title = pathway_name,
        x     = NULL,
        y     = "Enrichment Score"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(1, N),
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor  = ggplot2::element_blank(),
        panel.grid.major  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
        plot.title        = ggplot2::element_text(face = "bold", size = 11),
        axis.title.y      = ggplot2::element_text(face = "bold"),
        axis.text.x       = ggplot2::element_blank(),
        axis.ticks.x      = ggplot2::element_blank()
      )

    # --- Panel 2: gene hit rug ----------------------------------------------
    # geom_segment is used so that each mark renders at least 1px wide
    # regardless of how many genes are in the rank list (geom_tile with
    # width = 1 data-unit becomes sub-pixel for large N and disappears).
    hit_pos <- which(in_set)
    rug_df  <- data.frame(x = hit_pos)

    p_rug <- ggplot2::ggplot(rug_df, ggplot2::aes(x = x)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = x, y = 0, yend = 1),
        colour    = "black",
        linewidth = 0.3
      ) +
      ggplot2::scale_x_continuous(
        limits = c(1, N),
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = NULL,
        y = sprintf("Hits\n(n = %d)", k)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid        = ggplot2::element_blank(),
        axis.text         = ggplot2::element_blank(),
        axis.ticks        = ggplot2::element_blank(),
        axis.title.y      = ggplot2::element_text(face = "bold", size = 9)
      )

    # --- Panel 3: rank metric bar -------------------------------------------
    stat_df <- data.frame(x = seq_len(N), stat = stats)

    p_stat <- ggplot2::ggplot(stat_df, ggplot2::aes(x = x, y = stat)) +
      ggplot2::geom_col(
        ggplot2::aes(fill = stat > 0),
        width = 1,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC")) +
      ggplot2::scale_x_continuous(
        limits = c(1, N),
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::labs(
        x = "Gene rank",
        y = "Rank metric"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
        axis.title         = ggplot2::element_text(face = "bold"),
        axis.title.y       = ggplot2::element_text(face = "bold", size = 9)
      )

    # --- Compose with cowplot (align = "v", axis = "lr") -------------------
    # cowplot equalises the widths of left/right axis regions across all three
    # panels so the inner drawing areas are perfectly x-aligned.
    cowplot::plot_grid(
      p_es, p_rug, p_stat,
      ncol        = 1,
      align       = "v",
      axis        = "lr",
      rel_heights = c(0.55, 0.15, 0.30)
    )
  }

  # ---- output --------------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  plot_list <- lapply(names(pathways), function(pw_name) {
    p         <- .one_pathway_plot(pw_name, pathways[[pw_name]], stats)
    safe_name <- gsub("[^A-Za-z0-9._-]", "_", pw_name)
    out_file  <- file.path(out_dir, paste0(safe_name, ".png"))
    ggplot2::ggsave(
      filename = out_file,
      plot     = p,
      device   = "png",
      width    = width,
      height   = height,
      dpi      = 150
    )
    p
  })
  names(plot_list) <- names(pathways)

  invisible(plot_list)
}
