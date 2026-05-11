################################################################################
# modules/qc_visualization.R — Pre-filter QC visualizations
#
# Generates QC plots from the full unfiltered merged Seurat object.
# All grouping variables are derived automatically from METADATA_COVARIATE_COLS.
#
# Expects in environment:
#   seurat_obj, SAMPLE_METADATA, METADATA_COVARIATE_COLS, SPLIT_BY,
#   MITO_COL, MITO_PATTERN, QC_MAX_DECONTX, QC_MIN_FEATURES,
#   QC_MAX_FEATURES, QC_MAX_MITO_PCT,
#   stamp_pdf()
#
# Side effects:
#   Adds MITO_COL to seurat_obj metadata (in place)
#   Writes PDF plots to PLOT_DIR
################################################################################

run_qc_visualization <- function(seurat_obj) {

  message("\n── QC visualizations ────────────────────────────────────────────────────────")

  # Add mitochondrial % if not already present
  if (!MITO_COL %in% colnames(seurat_obj@meta.data))
    seurat_obj[[MITO_COL]] <- PercentageFeatureSet(seurat_obj, pattern = MITO_PATTERN)

  QC_METRICS <- c("nCount_RNA", "nFeature_RNA", MITO_COL, "decontX_contamination")
  QC_LABELS  <- c("UMI counts", "Genes detected", "Mitochondrial %", "DecontX contamination")
  GROUP_VARS <- c(SPLIT_BY, METADATA_COVARIATE_COLS)

  qc_df <- seurat_obj@meta.data[, c(QC_METRICS, GROUP_VARS), drop = FALSE]

  qc_theme <- function() {
    theme_bw(base_size = 11) +
      theme(axis.text.x    = element_text(angle = 45, hjust = 1),
            strip.text      = element_text(face = "bold"),
            legend.position = "none")
  }

  # ── 1. Violin + boxplot: each QC metric × each grouping variable ─────────────

  for (grp in GROUP_VARS) {
    n_groups <- length(unique(qc_df[[grp]]))
    pdf(stamp_pdf("QC_violin", grp), height = 5, width = max(6, n_groups * 0.8))
    for (k in seq_along(QC_METRICS)) {
      p <- ggplot(qc_df, aes(x = .data[[grp]], y = .data[[QC_METRICS[k]]],
                              fill = .data[[grp]])) +
        geom_violin(scale = "width", trim = TRUE) +
        geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.7) +
        labs(title = paste(QC_LABELS[k], "by", grp),
             x = grp, y = QC_LABELS[k]) +
        qc_theme()
      print(p)
    }
    dev.off()
    message("  Violin plots saved: QC_violin_", grp)
  }

  # ── 2. DecontX contamination density curves by grouping variable ──────────────

  for (grp in GROUP_VARS) {
    p <- ggplot(qc_df, aes(x = decontX_contamination,
                            color = .data[[grp]], fill = .data[[grp]])) +
      geom_density(alpha = 0.2) +
      geom_vline(xintercept = QC_MAX_DECONTX, linetype = "dashed",
                 color = "red", linewidth = 0.7) +
      annotate("text", x = QC_MAX_DECONTX + 0.01, y = Inf, vjust = 1.5,
               hjust = 0, size = 3, color = "red",
               label = paste("threshold =", QC_MAX_DECONTX)) +
      labs(title = paste("DecontX contamination density by", grp),
           x = "DecontX contamination score", y = "Density",
           color = grp, fill = grp) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    pdf(stamp_pdf("QC_decontX_density", grp), height = 4, width = 6)
    print(p)
    dev.off()
    message("  DecontX density saved: QC_decontX_density_", grp)
  }

  # ── 3. Pairwise scatter plots coloured by DecontX score ───────────────────────

  scatter_pairs <- list(
    list(x = "nCount_RNA",   y = "nFeature_RNA",         xl = "UMI counts",      yl = "Genes detected"),
    list(x = "nCount_RNA",   y = MITO_COL,               xl = "UMI counts",      yl = "Mitochondrial %"),
    list(x = "nCount_RNA",   y = "decontX_contamination", xl = "UMI counts",      yl = "DecontX contamination"),
    list(x = "nFeature_RNA", y = "decontX_contamination", xl = "Genes detected",  yl = "DecontX contamination"),
    list(x = MITO_COL,       y = "decontX_contamination", xl = "Mitochondrial %", yl = "DecontX contamination")
  )

  plot_idx      <- sample(nrow(qc_df), min(nrow(qc_df), 50000))
  scatter_plots <- lapply(scatter_pairs, function(pair) {
    p <- ggplot(qc_df[plot_idx, ],
                aes(x = .data[[pair$x]], y = .data[[pair$y]],
                    color = decontX_contamination)) +
      geom_point(size = 0.2, alpha = 0.4) +
      scale_color_viridis_c(option = "magma", name = "DecontX") +
      labs(x = pair$xl, y = pair$yl) +
      theme_bw(base_size = 10)
    if (pair$x == "nFeature_RNA")
      p <- p +
        geom_vline(xintercept = QC_MIN_FEATURES, linetype = "dashed", color = "blue",  linewidth = 0.5) +
        geom_vline(xintercept = QC_MAX_FEATURES, linetype = "dashed", color = "blue",  linewidth = 0.5)
    if (pair$y == "nFeature_RNA")
      p <- p +
        geom_hline(yintercept = QC_MIN_FEATURES, linetype = "dashed", color = "blue",  linewidth = 0.5) +
        geom_hline(yintercept = QC_MAX_FEATURES, linetype = "dashed", color = "blue",  linewidth = 0.5)
    if (pair$y == MITO_COL)
      p <- p + geom_hline(yintercept = QC_MAX_MITO_PCT,  linetype = "dashed", color = "red",   linewidth = 0.5)
    if (pair$y == "decontX_contamination")
      p <- p + geom_hline(yintercept = QC_MAX_DECONTX,   linetype = "dashed", color = "red",   linewidth = 0.5)
    p
  })

  pdf(stamp_pdf("QC_scatter_metrics"), height = 4, width = 4 * length(scatter_pairs))
  print(Reduce(`|`, scatter_plots))
  dev.off()
  message("  Scatter plots saved: QC_scatter_metrics")

  # ── 4. Cell count bar charts: per sample and per covariate ────────────────────

  pdf(stamp_pdf("QC_cell_counts"), height = 5,
      width = max(6, nrow(SAMPLE_METADATA) * 0.7))

  count_df          <- as.data.frame(table(qc_df[[SPLIT_BY]]))
  colnames(count_df) <- c("sample_id", "n_cells")
  print(
    ggplot(count_df, aes(x = sample_id, y = n_cells, fill = sample_id)) +
      geom_col() +
      geom_text(aes(label = n_cells), vjust = -0.3, size = 3) +
      labs(title = "Cells per sample (pre-filter)", x = NULL, y = "Cell count") +
      qc_theme()
  )

  for (grp in METADATA_COVARIATE_COLS) {
    grp_df            <- as.data.frame(table(qc_df[[grp]]))
    colnames(grp_df)  <- c("group", "n_cells")
    print(
      ggplot(grp_df, aes(x = group, y = n_cells, fill = group)) +
        geom_col() +
        geom_text(aes(label = n_cells), vjust = -0.3, size = 3) +
        labs(title = paste("Cells per", grp, "(pre-filter)"), x = grp, y = "Cell count") +
        qc_theme()
    )
  }
  dev.off()
  message("  Cell count bar charts saved: QC_cell_counts")

  # ── 5. QC metric correlation heatmap ──────────────────────────────────────────

  cor_mat           <- cor(qc_df[, QC_METRICS], use = "pairwise.complete.obs")
  cor_long          <- reshape2::melt(cor_mat)
  colnames(cor_long) <- c("Var1", "Var2", "Correlation")
  label_map         <- setNames(QC_LABELS, QC_METRICS)
  cor_long$Var1     <- label_map[as.character(cor_long$Var1)]
  cor_long$Var2     <- label_map[as.character(cor_long$Var2)]

  pdf(stamp_pdf("QC_metric_correlation_heatmap"), height = 4, width = 5)
  print(
    ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(Correlation, 2)), size = 3.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(title = "QC metric correlations (pre-filter)", x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  message("  Correlation heatmap saved: QC_metric_correlation_heatmap")

  message("QC visualizations complete.\n")
  invisible(seurat_obj)
}
