################################################################################
# modules/seurat_qc.R — Seurat QC filtering and clustering
#
# Applies QC thresholds to singlets, normalizes with SCTransform,
# runs PCA/UMAP/clustering, and saves the final Seurat object.
#
# Expects in environment:
#   MITO_COL, MITO_PATTERN, SPLIT_BY, PROJECT_NAME,
#   QC_MIN_FEATURES, QC_MAX_FEATURES, QC_MAX_MITO_PCT, QC_MAX_DECONTX,
#   QC_MAX_MITO_PCT, QC_MAX_DECONTX,
#   SCT_VARS_TO_REGRESS, N_PCS,
#   drop_zero_umi(), select_pcs(), stamp(), stamp_pdf()
#
# Returns:
#   seurat_filt — final filtered, normalized, and clustered Seurat object
################################################################################

run_seurat_qc <- function(seurat_singlets) {

  message("\n── Seurat QC filtering ──────────────────────────────────────────────────────")

  # Ensure mito % is present (should carry through from earlier steps)
  if (!MITO_COL %in% colnames(seurat_singlets@meta.data))
    seurat_singlets[[MITO_COL]] <- PercentageFeatureSet(seurat_singlets,
                                                         pattern = MITO_PATTERN)

  # ── Pre-filter QC plots ────────────────────────────────────────────────────────

  pdf(stamp_pdf("QC_violinPlots_PRE_FILTER"), height = 5, width = 10)
  print(VlnPlot(seurat_singlets,
                features = c("nFeature_RNA", "nCount_RNA", MITO_COL, "decontX_contamination"),
                ncol = 2, pt.size = 0, group.by = SPLIT_BY))
  dev.off()

  plot1 <- FeatureScatter(seurat_singlets, feature1 = "nCount_RNA", feature2 = MITO_COL)
  plot2 <- FeatureScatter(seurat_singlets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(seurat_singlets, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

  pdf(stamp_pdf("QC_scatter_PRE_FILTER"), height = 5, width = 12)
  print(plot1 | plot2 | plot3)
  dev.off()
  message("  Pre-filter QC plots saved.")

  # ── Apply QC thresholds ────────────────────────────────────────────────────────
  # Built as a logical index rather than subset() to support dynamic MITO_COL name

  # Default QC_MAX_FEATURES for older configs that predate this variable
  if (!exists("QC_MAX_FEATURES")) QC_MAX_FEATURES <- Inf

  md   <- seurat_singlets@meta.data
  keep <- md$nFeature_RNA          > QC_MIN_FEATURES &
          md$nFeature_RNA          < QC_MAX_FEATURES &
          md[[MITO_COL]]           < QC_MAX_MITO_PCT  &
          md$decontX_contamination < QC_MAX_DECONTX

  seurat_filt <- seurat_singlets[, keep]

  message("  Cells after QC filtering: ", ncol(seurat_filt))
  message("  Cell counts per sample:")
  print(table(seurat_filt@meta.data[[SPLIT_BY]]))

  # ── Post-filter QC plots ───────────────────────────────────────────────────────

  pdf(stamp_pdf("QC_violinPlots_POST_FILTER"), height = 5, width = 10)
  print(VlnPlot(seurat_filt,
                features = c("nFeature_RNA", "nCount_RNA", MITO_COL, "decontX_contamination"),
                ncol = 2, pt.size = 0, group.by = SPLIT_BY))
  dev.off()
  message("  Post-filter QC plots saved.")

  save(seurat_filt, file = stamp("Seurat_post_QC", ext = "RData"))

  # ── Normalize, SCTransform, PCA ───────────────────────────────────────────────

  seurat_filt <- NormalizeData(seurat_filt, normalization.method = "LogNormalize",
                                scale.factor = 10000)
  seurat_filt <- drop_zero_umi(seurat_filt, "final")
  seurat_filt <- SCTransform(seurat_filt,
                              vars.to.regress = filter_regression_vars(seurat_filt, SCT_VARS_TO_REGRESS),
                              verbose = FALSE)
  seurat_filt <- RunPCA(seurat_filt, npcs = N_PCS, verbose = FALSE)

  # ── PC selection and elbow plot ───────────────────────────────────────────────

  message("  PC selection for final clustering:")
  pcs_final <- select_pcs(seurat_filt)

  pct     <- seurat_filt[["pca"]]@stdev / sum(seurat_filt[["pca"]]@stdev) * 100
  cumu    <- cumsum(pct)
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))

  pdf(stamp_pdf("ElbowPlot_threshold"), height = 5, width = 6)
  print(
    ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs_final)) +
      geom_text() +
      geom_vline(xintercept = 90, color = "grey") +
      geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      theme_bw() +
      labs(title = "PC selection", x = "Cumulative % variance", y = "% variance per PC")
  )
  dev.off()

  # ── UMAP and clustering ───────────────────────────────────────────────────────

  seurat_filt <- RunUMAP(seurat_filt, dims = 1:pcs_final,
                          umap.method = UMAP_METHOD, metric = UMAP_METRIC, verbose = FALSE)
  seurat_filt <- FindNeighbors(seurat_filt, dims = 1:pcs_final, verbose = FALSE)

  # ── Resolution selection ──────────────────────────────────────────────────────
  # If CLUSTER_RESOLUTION is set: use it directly.
  # If NULL: sweep a range of resolutions with clustree and pick automatically.

  # Defaults for older configs that don't define these variables
  if (!exists("CLUSTER_RESOLUTION"))  CLUSTER_RESOLUTION  <- NULL
  if (!exists("CLUSTREE_RESOLUTIONS")) CLUSTREE_RESOLUTIONS <- seq(0.2, 2.0, by = 0.2)

  if (!is.null(CLUSTER_RESOLUTION)) {

    # ── Fixed resolution ────────────────────────────────────────────────────────
    message("  Using fixed clustering resolution: ", CLUSTER_RESOLUTION)
    seurat_filt <- FindClusters(seurat_filt, resolution = CLUSTER_RESOLUTION,
                                 verbose = FALSE)
    chosen_res  <- CLUSTER_RESOLUTION

  } else {

    # ── Clustree auto-selection ─────────────────────────────────────────────────
    message("  Running clustree resolution sweep: ",
            paste(CLUSTREE_RESOLUTIONS, collapse = ", "))

    # Run FindClusters at every resolution in the sweep
    seurat_filt <- FindClusters(seurat_filt, resolution = CLUSTREE_RESOLUTIONS,
                                 verbose = FALSE)

    # Build clustree and save the plot
    pdf(stamp_pdf("clustree_resolution_sweep"), height = 10, width = 8)
    print(clustree(seurat_filt, prefix = "SCT_snn_res."))
    dev.off()
    message("  Clustree plot saved: clustree_resolution_sweep")

    # Auto-select the resolution where the number of new clusters first
    # stabilises — defined as the lowest resolution where no further splits occur
    # (i.e. the first resolution after which cluster count stops increasing).
    res_col_prefix <- "SCT_snn_res."
    res_cols       <- grep(res_col_prefix, colnames(seurat_filt@meta.data), value = TRUE)
    n_clusters     <- sapply(res_cols, function(col)
                               length(unique(seurat_filt@meta.data[[col]])))

    # Find the first resolution where cluster count stops increasing
    delta        <- diff(n_clusters)
    stable_idx   <- which(delta == 0)[1]
    chosen_res   <- if (!is.na(stable_idx)) {
                     as.numeric(sub(res_col_prefix, "", res_cols[stable_idx]))
                   } else {
                     # If cluster count never fully stabilises, pick the elbow:
                     # resolution where the rate of increase drops the most
                     elbow_idx  <- which.min(diff(delta)) + 1
                     as.numeric(sub(res_col_prefix, "", res_cols[elbow_idx]))
                   }

    message("  Clustree selected resolution: ", chosen_res)

    # Set the chosen resolution as the active cluster identity
    Idents(seurat_filt) <- paste0(res_col_prefix, chosen_res)
    seurat_filt@meta.data$seurat_clusters <- seurat_filt@meta.data[[paste0(res_col_prefix, chosen_res)]]

  }

  message("  Number of clusters: ",
          length(unique(seurat_filt@meta.data$seurat_clusters)))

  pdf(stamp_pdf(paste0("UMAP_SeuratClusters_res", chosen_res)), height = 5, width = 7)
  print(DimPlot(seurat_filt, reduction = "umap", label = TRUE) + NoLegend() +
          ggtitle(paste("Resolution:", chosen_res)))
  dev.off()

  save(seurat_filt, file = stamp("Seurat_final", ext = "RData"))
  message("  Final Seurat object saved.")

  invisible(seurat_filt)
}
