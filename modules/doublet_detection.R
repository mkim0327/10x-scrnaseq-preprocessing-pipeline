################################################################################
# modules/doublet_detection.R — Doublet detection
#
# Runs DoubletFinder and scds (cxds_bcds_hybrid) on each sample separately,
# combines calls using a union rule, visualizes results, and returns the
# merged Seurat object with singlets only.
#
# Doublets are detected BEFORE QC filtering so that the full cell neighborhood
# context is available to both callers.
#
# Expects in environment:
#   seurat_obj, SAMPLE_METADATA, SPLIT_BY, PROJECT_NAME,
#   SCT_VARS_TO_REGRESS, N_PCS, DOUBLET_EXTRA_RATE, QC_MIN_CELLS,
#   MITO_COL, MITO_PATTERN, N_WORKERS,
#   DOUBLET_METHOD — one of: "union" (default), "intersection",
#                            "DoubletFinder", "scds"
#   make_doublet_rate_model(), drop_zero_umi(), select_pcs(),
#   stamp(), stamp_pdf()
#
# Returns:
#   seurat_singlets — Seurat object containing singlets only
################################################################################

run_doublet_detection <- function(seurat_obj) {

  message("\n── Doublet detection ────────────────────────────────────────────────────────")

  # ── 3a. Minimal preprocessing required by doublet callers ────────────────────

  # Remove genes expressed in very few cells
  num_cells      <- Matrix::rowSums(seurat_obj@assays$RNA@counts > 0)
  genes_use      <- names(num_cells[num_cells >= QC_MIN_CELLS])
  seurat_prefilt <- subset(seurat_obj, features = genes_use)

  if (!MITO_COL %in% colnames(seurat_prefilt@meta.data))
    seurat_prefilt[[MITO_COL]] <- PercentageFeatureSet(seurat_prefilt, pattern = MITO_PATTERN)

  seurat_prefilt <- drop_zero_umi(seurat_prefilt, "prefilt")
  seurat_prefilt <- SCTransform(seurat_prefilt,
                                 vars.to.regress = filter_regression_vars(seurat_prefilt, SCT_VARS_TO_REGRESS),
                                 verbose = FALSE)
  seurat_prefilt <- RunPCA(seurat_prefilt, npcs = N_PCS, verbose = FALSE)

  message("PC selection for doublet detection:")
  pcs <- select_pcs(seurat_prefilt)

  seurat_prefilt <- RunUMAP(seurat_prefilt, dims = 1:pcs,
                             umap.method = UMAP_METHOD, metric = UMAP_METRIC, verbose = FALSE)
  seurat_prefilt <- FindNeighbors(seurat_prefilt, dims = 1:pcs, verbose = FALSE)
  # Resolution here is intentionally fixed — this clustering is only used to
  # estimate homotypic doublet proportions for DoubletFinder, not for biology.
  # The final biological clustering resolution is set in config.R.
  seurat_prefilt <- FindClusters(seurat_prefilt, resolution = 2, verbose = FALSE)

  # ── 3b. DoubletFinder ────────────────────────────────────────────────────────

  pred_dblt_lm   <- make_doublet_rate_model()
  sample_list    <- SplitObject(seurat_prefilt, split.by = SPLIT_BY)

  pred_dblt_rate <- predict(
    pred_dblt_lm,
    data.frame(cell_number = 1.1 * sapply(sample_list, ncol))
  ) / 100

  message("  Running DoubletFinder on ", length(sample_list),
          " sample(s) with ", N_WORKERS, " worker(s)...")

  df_results <- future.apply::future_lapply(
    seq_along(sample_list),
    function(i) {
      obj <- sample_list[[i]]
      s   <- names(sample_list)[i]

      # Per-sample SCTransform: batch covariates are constant within a sample
      sample_vars <- filter_regression_vars(obj, SCT_VARS_TO_REGRESS)
      obj <- SCTransform(obj, vars.to.regress = sample_vars, verbose = FALSE)

      # num.cores = 1 avoids nested parallelism with the outer future workers
      sweep_res   <- paramSweep_v3(obj, PCs = 1:pcs, sct = TRUE, num.cores = 1)
      sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
      bcmvn       <- find.pK(sweep_stats)
      pk          <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

      nExp <- round((pred_dblt_rate[i] + DOUBLET_EXTRA_RATE) * ncol(obj))
      obj  <- doubletFinder_v3(obj, PCs = 1:pcs, pN = 0.25, pK = pk,
                                nExp = nExp, reuse.pANN = FALSE, sct = TRUE)

      # Standardise column name
      df_col <- grep("DF.classifications_0.25", colnames(obj@meta.data), value = TRUE)
      colnames(obj@meta.data)[colnames(obj@meta.data) == df_col] <- "DoubletFinder"
      obj
    },
    future.seed = TRUE
  )
  sample_list <- setNames(df_results, names(sample_list))

  pdf(stamp_pdf("DoubletFinder_UMAPs"), height = 5, width = 6)
  for (i in seq_along(sample_list)) {
    print(DimPlot(sample_list[[i]], reduction = "umap", group.by = "DoubletFinder") +
            ggtitle(names(sample_list)[i]))
  }
  dev.off()
  message("  DoubletFinder UMAPs saved.")

  # ── 3c. scds hybrid scoring ───────────────────────────────────────────────────
  # Convert to SCE using RNA assay (not SCT) so scds operates on raw counts

  message("  Running scds on ", length(sample_list),
          " sample(s) with ", N_WORKERS, " worker(s)...")

  scds_results <- future.apply::future_lapply(
    seq_along(sample_list),
    function(i) {
      obj                    <- sample_list[[i]]
      DefaultAssay(obj)      <- "RNA"
      sce                    <- as.SingleCellExperiment(obj)
      sce                    <- cxds_bcds_hybrid(sce)
      n_db                   <- round(pred_dblt_rate[i] * ncol(sce))
      top_db                 <- order(sce$hybrid_score, decreasing = TRUE)[seq_len(n_db)]
      sce$scds               <- "Singlet"
      sce$scds[top_db]       <- "Doublet"
      sce
    },
    future.seed = TRUE
  )
  scds_list <- setNames(scds_results, names(sample_list))

  # Plot scds scores on the Seurat UMAP (from the DoubletFinder sample_list)
  # rather than trying to use a UMAP from the SCE (which doesn't have one)
  pdf(stamp_pdf("scds_UMAPs"), height = 5, width = 6)
  for (i in seq_along(sample_list)) {
    s          <- names(sample_list)[i]
    scds_calls <- scds_list[[i]]$scds
    names(scds_calls) <- colnames(scds_list[[i]])
    sample_list[[i]]@meta.data$scds_plot <- scds_calls[colnames(sample_list[[i]])]
    print(DimPlot(sample_list[[i]], reduction = "umap", group.by = "scds_plot") +
            ggtitle(s) + labs(color = "scds"))
  }
  dev.off()
  message("  scds UMAPs saved.")

  # ── 3d. Merge per-sample results and transfer scds scores ─────────────────────

  if (length(sample_list) == 1) {
    seurat_doublets <- sample_list[[1]]
  } else {
    seurat_doublets <- merge(sample_list[[1]], y = sample_list[-1], project = PROJECT_NAME)
  }

  # Remove per-lane pANN columns (DoubletFinder artefact)
  pann_cols <- grep("pANN", colnames(seurat_doublets@meta.data))
  if (length(pann_cols) > 0)
    seurat_doublets@meta.data <- seurat_doublets@meta.data[, -pann_cols]

  # Transfer scds labels using colData indexing ($ on SCE doesn't support
  # character vector subsetting for cell names)
  seurat_doublets@meta.data$scds_hybrid <- NA_character_
  for (i in seq_along(scds_list)) {
    shared <- intersect(colnames(scds_list[[i]]), colnames(seurat_doublets))
    seurat_doublets@meta.data[shared, "scds_hybrid"] <-
      as.character(colData(scds_list[[i]])[shared, "scds"])
  }

  # ── 3e. Combine doublet calls per selected method ────────────────────────────

  valid_methods <- c("union", "intersection", "DoubletFinder", "scds")
  if (!exists("DOUBLET_METHOD") || !DOUBLET_METHOD %in% valid_methods) {
    warning("DOUBLET_METHOD not set or invalid. Defaulting to 'union'. ",
            "Valid options: ", paste(valid_methods, collapse = ", "))
    DOUBLET_METHOD <- "union"
  }

  seurat_doublets@meta.data$DoubletCall <- switch(
    DOUBLET_METHOD,
    "union" = ifelse(
      bitOr(seurat_doublets@meta.data$DoubletFinder == "Doublet",
            seurat_doublets@meta.data$scds_hybrid   == "Doublet") > 0,
      "Doublet", "Singlet"
    ),
    "intersection" = ifelse(
      seurat_doublets@meta.data$DoubletFinder == "Doublet" &
      seurat_doublets@meta.data$scds_hybrid   == "Doublet",
      "Doublet", "Singlet"
    ),
    "DoubletFinder" = seurat_doublets@meta.data$DoubletFinder,
    "scds"          = seurat_doublets@meta.data$scds_hybrid
  )

  message("  Doublet method: ", DOUBLET_METHOD)
  message("  Doublet call summary:")
  print(table(seurat_doublets@meta.data$DoubletCall))
  message("  DoubletFinder vs scds cross-tabulation:")
  print(table(DoubletFinder = seurat_doublets@meta.data$DoubletFinder,
              scds          = seurat_doublets@meta.data$scds_hybrid))

  # Visualise doublet calls on UMAP
  seurat_doublets <- drop_zero_umi(seurat_doublets, "doublets")
  seurat_doublets <- SCTransform(seurat_doublets,
                                  vars.to.regress = filter_regression_vars(seurat_doublets, SCT_VARS_TO_REGRESS),
                                  verbose = FALSE)
  seurat_doublets <- RunPCA(seurat_doublets,  npcs = N_PCS,  verbose = FALSE)
  seurat_doublets <- RunUMAP(seurat_doublets, dims = 1:pcs,
                              umap.method = UMAP_METHOD, metric = UMAP_METRIC, verbose = FALSE)

  pdf(stamp_pdf("DoubletCall_UMAP"), height = 5, width = 6)
  print(DimPlot(seurat_doublets, reduction = "umap", group.by = "DoubletCall"))
  dev.off()

  save(seurat_doublets, file = stamp("Seurat_AnnotatedDoublets", ext = "RData"))

  # ── 3f. Retain singlets only ──────────────────────────────────────────────────

  seurat_singlets <- subset(seurat_doublets, subset = DoubletCall == "Singlet")
  message("  Cells after doublet removal: ", ncol(seurat_singlets))

  invisible(seurat_singlets)
}
