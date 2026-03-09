plot_pca <- function(expr_mat, metadata, color_var, out_dir = "figures/pca", width = 7, height = 5) {
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
    stop("`metadata` must contain a `Sample` column matching expression matrix sample IDs.")
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
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (is.numeric(color_data)) {
    plot_obj <- plot_obj +
      ggplot2::scale_color_gradientn(colors = grDevices::hcl.colors(9, "TealGrn"))
  } else {
    n_levels <- length(unique(as.character(color_data)))
    plot_obj <- plot_obj +
      ggplot2::scale_color_manual(values = grDevices::hcl.colors(n_levels, "Dark 3"))
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", color_var)
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

plot_density_counts <- function(expr_mat, metadata, color_var, out_dir = "figures/density", width = 8, height = 6, is_raw_counts = NULL) {
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

  if (!is.null(is_raw_counts) && (!is.logical(is_raw_counts) || length(is_raw_counts) != 1)) {
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
    stop("`metadata` must contain a `Sample` column matching count matrix sample IDs.")
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
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (is.numeric(color_data)) {
    plot_obj <- plot_obj +
      ggplot2::scale_color_gradientn(colors = grDevices::hcl.colors(9, "TealGrn"))
  } else {
    n_levels <- length(unique(as.character(color_data)))
    plot_obj <- plot_obj +
      ggplot2::scale_color_manual(values = grDevices::hcl.colors(n_levels, "Dark 3"))
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", color_var)
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

plot_ma <- function(counts_mat, out_dir = "figures/ma", file_name = "ma_plots.pdf", plots_per_page = 10, ncol = 2, width = 11, height = 14) {
  if (!is.matrix(counts_mat)) {
    stop("`counts_mat` must be a matrix.")
  }

  if (!is.numeric(counts_mat)) {
    stop("`counts_mat` must be a numeric matrix.")
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

  if (!is.numeric(plots_per_page) || length(plots_per_page) != 1 || plots_per_page <= 0) {
    stop("`plots_per_page` must be a single positive number.")
  }

  if (!is.numeric(ncol) || length(ncol) != 1 || ncol <= 0) {
    stop("`ncol` must be a single positive number.")
  }

  sample_ids <- colnames(counts_mat)
  if (is.null(sample_ids)) {
    stop("`counts_mat` must have column names corresponding to sample IDs.")
  }

  log_mat <- log2(counts_mat + 1)
  ref_log <- rowMeans(log_mat, na.rm = TRUE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  plot_list <- vector("list", length(sample_ids))
  names(plot_list) <- sample_ids

  for (i in seq_along(sample_ids)) {
    sample_name <- sample_ids[i]
    sample_log <- log_mat[, i]

    plot_df <- data.frame(
      A = 0.5 * (sample_log + ref_log),
      M = sample_log - ref_log,
      stringsAsFactors = FALSE
    )
    plot_df <- plot_df[is.finite(plot_df$A) & is.finite(plot_df$M), , drop = FALSE]

    plot_obj <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = A, y = M)
    ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "grey85") +
      ggplot2::geom_point(size = 1.1, alpha = 0.55, color = "#0072B2") +
      ggplot2::labs(
        title = "MA Plot",
        subtitle = sprintf("Sample: %s", sample_name),
        x = "A: mean log2 expression",
        y = "M: log2 fold-change vs average sample"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      )

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
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow, ncol = ncol)))

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
