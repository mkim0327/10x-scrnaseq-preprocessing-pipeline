################################################################################
# modules/cellranger_summary.R — CellRanger experiment summary
#
# Reads metrics_summary.csv from each sample's CellRanger output, saves a
# merged summary table, and generates visualizations.
#
# Expects in environment:
#   SAMPLE_METADATA, METADATA_COVARIATE_COLS, CELLRANGER_OUTDIR,
#   stamp(), stamp_pdf()
#
# Returns:
#   qc_stats  — data frame of raw CellRanger metrics (one row per sample)
################################################################################

run_cellranger_summary <- function() {

  message("\n── CellRanger experiment summary ────────────────────────────────────────────")

  # ── Read metrics files ───────────────────────────────────────────────────────

  metrics_files <- setNames(
    sapply(SAMPLE_METADATA$sample_id, function(s) cellranger_paths(s)$metrics),
    SAMPLE_METADATA$sample_id
  )

  missing_metrics <- metrics_files[!file.exists(metrics_files)]
  if (length(missing_metrics) > 0) {
    stop("metrics_summary.csv not found for:\n",
         paste(names(missing_metrics), collapse = "\n"),
         "\nHas CellRanger finished for all samples?")
  }

  qc_stats <- data.frame(
    do.call(rbind, lapply(metrics_files, read.csv)),
    row.names = SAMPLE_METADATA$sample_id
  )

  write.table(qc_stats, file = stamp("CellRanger_QC_summary", ext = "txt"),
              sep = "\t", row.names = TRUE, quote = FALSE)
  message("CellRanger QC summary saved.")
  print(qc_stats)

  # ── Prepare plotting data frame ──────────────────────────────────────────────
  # Convert % strings to numeric; strip thousands separators

  cr_plot_df <- as.data.frame(lapply(qc_stats, function(col) {
    cleaned <- gsub(",", "", as.character(col))
    cleaned <- gsub("%", "", cleaned)
    suppressWarnings(as.numeric(cleaned))
  }))
  cr_plot_df$sample_id <- rownames(qc_stats)

  # Join metadata covariates if any exist — guard against empty covariate list
  if (length(METADATA_COVARIATE_COLS) > 0) {
    cr_plot_df <- merge(
      cr_plot_df,
      SAMPLE_METADATA[, c("sample_id", METADATA_COVARIATE_COLS), drop = FALSE],
      by = "sample_id", all.x = TRUE
    )
  }

  # Classify metrics: percentage vs count
  pct_metrics <- names(qc_stats)[sapply(qc_stats, function(col)
    any(grepl("%", as.character(col))))]

  # Clean display names (dots from make.names → spaces)
  clean_name <- function(x) gsub("\\.", " ", x)

  message("\n── CellRanger summary visualizations ────────────────────────────────────────")

  # ── Bar charts: each metric × each grouping variable ────────────────────────

  for (grp in c("sample_id", METADATA_COVARIATE_COLS)) {
    n_groups <- length(unique(cr_plot_df[[grp]]))
    pdf(stamp_pdf("CellRanger_summary_by", grp),
        height = 4, width = max(5, n_groups * 0.9))
    for (metric in colnames(qc_stats)) {
      y_label <- if (metric %in% pct_metrics)
                   paste0(clean_name(metric), " (%)")
                 else
                   clean_name(metric)
      p <- ggplot(cr_plot_df,
                  aes(x = .data[["sample_id"]], y = .data[[metric]],
                      fill = .data[[grp]])) +
        geom_col(position = "dodge") +
        geom_text(aes(label = round(.data[[metric]], 1)), vjust = -0.3, size = 2.8) +
        labs(title = clean_name(metric), x = NULL, y = y_label, fill = grp) +
        theme_bw(base_size = 11) +
        theme(axis.text.x    = element_text(angle = 45, hjust = 1),
              legend.position = if (grp == "sample_id") "none" else "right")
      print(p)
    }
    dev.off()
    message("  Bar charts saved: CellRanger_summary_by_", grp)
  }

  # ── Z-scored heatmap: all metrics × all samples ──────────────────────────────

  numeric_cols <- intersect(colnames(cr_plot_df), colnames(qc_stats))
  cr_mat       <- as.matrix(cr_plot_df[, numeric_cols])
  rownames(cr_mat) <- cr_plot_df$sample_id

  cr_mat_z           <- scale(cr_mat)
  cr_mat_z[is.nan(cr_mat_z)] <- 0

  cr_long         <- reshape2::melt(cr_mat_z)
  colnames(cr_long) <- c("sample_id", "metric", "z_score")
  cr_long$metric  <- clean_name(as.character(cr_long$metric))

  pdf(stamp_pdf("CellRanger_summary_heatmap"),
      height = max(4, ncol(cr_mat) * 0.35),
      width  = max(5, nrow(cr_mat) * 0.6))
  print(
    ggplot(cr_long, aes(x = sample_id, y = metric, fill = z_score)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(z_score, 1)), size = 2.8) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0, name = "Z-score") +
      labs(title = "CellRanger metrics (z-scored across samples)",
           x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  message("  Heatmap saved: CellRanger_summary_heatmap")

  # ── Sequencing quality overview: all % metrics side-by-side ─────────────────

  if (length(pct_metrics) > 0) {
    pct_long        <- reshape2::melt(cr_plot_df[, c("sample_id", pct_metrics)],
                                      id.vars      = "sample_id",
                                      variable.name = "metric",
                                      value.name    = "pct")
    pct_long$metric <- clean_name(as.character(pct_long$metric))

    pdf(stamp_pdf("CellRanger_sequencing_quality"), height = 5, width = 8)
    print(
      ggplot(pct_long, aes(x = sample_id, y = pct, fill = metric)) +
        geom_col(position = "dodge") +
        labs(title = "Sequencing quality metrics by sample",
             x = NULL, y = "Percentage (%)", fill = "Metric") +
        theme_bw(base_size = 11) +
        theme(axis.text.x    = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2))
    )
    dev.off()
    message("  Sequencing quality chart saved: CellRanger_sequencing_quality")
  }

  message("CellRanger summary complete.\n")
  invisible(qc_stats)
}
