################################################################################
# modules/batch_evaluation.R — Batch effect evaluation
#
# Runs after final clustering to assess whether batch effects remain in the
# data. Produces visual and quantitative outputs to inform the decision of
# whether integration is needed.
#
# Metrics:
#   - UMAP coloured by each batch variable (visual)
#   - Silhouette score per batch variable in PCA space (quantitative)
#     Score range: -1 to 1
#       ~1   = cells cluster tightly with their own batch → strong batch effect
#       ~0   = batches overlap → no batch effect / well mixed
#       < 0  = cells closer to another batch → over-correction or noise
#   - Per-cluster batch composition bar chart (visual)
#     Uneven composition across clusters flags batch-driven clustering.
#
# Expects in environment:
#   seurat_filt, BATCH_COLS, SPLIT_BY, PLOT_DIR,
#   stamp(), stamp_pdf()
#
# Returns: invisible(NULL) — results are written to disk only
################################################################################

run_batch_evaluation <- function(seurat_filt) {

  message("\n── Batch effect evaluation ──────────────────────────────────────────────────")

  if (!requireNamespace("cluster", quietly = TRUE))
    stop("Package 'cluster' is required for silhouette scoring. ",
         "Install with: install.packages('cluster')")

  # ── Collect PCA embeddings and metadata ───────────────────────────────────────

  pca_embed <- Embeddings(seurat_filt, reduction = "pca")
  meta      <- seurat_filt@meta.data

  # Results table — one row per batch variable
  silhouette_results <- data.frame(
    batch_variable   = character(),
    mean_silhouette  = numeric(),
    interpretation   = character(),
    stringsAsFactors = FALSE
  )

  for (batch_var in BATCH_COLS) {

    if (!batch_var %in% colnames(meta)) {
      message("  Skipping '", batch_var, "' — not found in metadata.")
      next
    }

    batch_labels <- meta[[batch_var]]

    # Silhouette requires at least 2 groups with >= 1 cell each
    if (length(unique(batch_labels)) < 2) {
      message("  Skipping '", batch_var, "' — only one unique value.")
      next
    }

    message("  Evaluating batch variable: ", batch_var)

    # ── UMAP coloured by batch variable ─────────────────────────────────────────

    p_umap <- DimPlot(seurat_filt, reduction = "umap",
                      group.by = batch_var, label = FALSE) +
      labs(title = paste("UMAP coloured by", batch_var),
           color = batch_var) +
      theme_bw(base_size = 11)

    pdf(stamp_pdf("batch_UMAP", batch_var), height = 5, width = 7)
    print(p_umap)
    dev.off()
    message("    UMAP saved: batch_UMAP_", batch_var)

    # ── Silhouette score in PCA space ────────────────────────────────────────────
    # Use a distance matrix on the first N PCs (capped at 50 for memory)
    # Subsample if > 20,000 cells to keep compute tractable

    max_cells <- 20000
    if (nrow(pca_embed) > max_cells) {
      message("    Subsampling to ", max_cells, " cells for silhouette computation.")
      idx         <- sample(nrow(pca_embed), max_cells)
      pca_sub     <- pca_embed[idx, seq_len(min(50, ncol(pca_embed))), drop = FALSE]
      labels_sub  <- batch_labels[idx]
    } else {
      pca_sub    <- pca_embed[, seq_len(min(50, ncol(pca_embed))), drop = FALSE]
      labels_sub <- batch_labels
    }

    dist_mat  <- dist(pca_sub)
    label_int <- as.integer(factor(labels_sub))
    sil       <- cluster::silhouette(label_int, dist_mat)
    sil_mean  <- round(mean(sil[, "sil_width"]), 4)

    interpretation <- dplyr::case_when(
      sil_mean >  0.25 ~ "Strong batch effect — integration recommended",
      sil_mean >  0.10 ~ "Moderate batch effect — integration may help",
      sil_mean >= 0.00 ~ "Weak batch effect — integration likely not needed",
      TRUE             ~ "Negative score — batches well mixed or noise"
    )

    message("    Mean silhouette score (", batch_var, "): ", sil_mean,
            " [", interpretation, "]")

    silhouette_results <- rbind(silhouette_results, data.frame(
      batch_variable  = batch_var,
      mean_silhouette = sil_mean,
      interpretation  = interpretation,
      stringsAsFactors = FALSE
    ))

    # ── Silhouette distribution plot ─────────────────────────────────────────────

    sil_df        <- as.data.frame(sil[, c("cluster", "sil_width")])
    sil_df$batch  <- labels_sub
    sil_df        <- sil_df[order(sil_df$cluster, sil_df$sil_width), ]
    sil_df$cell   <- seq_len(nrow(sil_df))

    p_sil <- ggplot(sil_df, aes(x = cell, y = sil_width, fill = batch)) +
      geom_col(width = 1) +
      geom_hline(yintercept = sil_mean, linetype = "dashed",
                 color = "black", linewidth = 0.6) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
      annotate("text", x = nrow(sil_df) * 0.02, y = sil_mean + 0.03,
               label = paste("mean =", sil_mean), hjust = 0, size = 3) +
      labs(title   = paste("Silhouette scores by", batch_var),
           x = "Cells (sorted by batch and score)",
           y = "Silhouette width",
           fill = batch_var) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    pdf(stamp_pdf("batch_silhouette", batch_var), height = 4, width = 7)
    print(p_sil)
    dev.off()
    message("    Silhouette plot saved: batch_silhouette_", batch_var)

    # ── Per-cluster batch composition bar chart ───────────────────────────────────
    # Uneven batch composition within clusters flags batch-driven clustering

    comp_df <- as.data.frame(table(
      Cluster = seurat_filt@meta.data$seurat_clusters,
      Batch   = meta[[batch_var]]
    ))
    # Normalise to proportion within each cluster
    comp_df <- comp_df |>
      dplyr::group_by(Cluster) |>
      dplyr::mutate(proportion = Freq / sum(Freq)) |>
      dplyr::ungroup()

    n_clusters <- length(unique(comp_df$Cluster))

    p_comp <- ggplot(comp_df, aes(x = Cluster, y = proportion, fill = Batch)) +
      geom_col(position = "stack") +
      geom_hline(yintercept = 1 / length(unique(comp_df$Batch)),
                 linetype = "dashed", color = "black", linewidth = 0.5) +
      labs(title = paste("Cluster composition by", batch_var),
           subtitle = "Dashed line = expected proportion if perfectly mixed",
           x = "Cluster", y = "Proportion", fill = batch_var) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(
              angle = if (n_clusters > 20) 45 else 0, hjust = 1))

    pdf(stamp_pdf("batch_cluster_composition", batch_var),
        height = 5, width = max(6, n_clusters * 0.4))
    print(p_comp)
    dev.off()
    message("    Cluster composition plot saved: batch_cluster_composition_", batch_var)
  }

  # ── Save summary table ────────────────────────────────────────────────────────

  if (nrow(silhouette_results) > 0) {
    write.table(silhouette_results,
                file = stamp("batch_evaluation_summary", ext = "txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    message("\n  ── Batch evaluation summary ─────────────────────────────────────")
    for (i in seq_len(nrow(silhouette_results))) {
      message(sprintf("  %-20s  silhouette = %6.4f  →  %s",
                      silhouette_results$batch_variable[i],
                      silhouette_results$mean_silhouette[i],
                      silhouette_results$interpretation[i]))
    }
    message("  ─────────────────────────────────────────────────────────────────")
    message("  Full results saved: batch_evaluation_summary.txt")
  }

  message("Batch evaluation complete.\n")
  invisible(NULL)
}
