################################################################################
# modules/decontx.R — Ambient RNA removal with DecontX
#
# Runs DecontX on each sample using the raw matrix as background, saves
# per-sample SCE objects, builds Seurat objects from corrected counts,
# drops zero-UMI cells introduced by rounding, merges all samples, and
# attaches metadata from SAMPLE_METADATA.
#
# Expects in environment:
#   SAMPLE_METADATA, METADATA_COVARIATE_COLS, PROJECT_NAME, N_WORKERS,
#   cellranger_paths(), stamp(), stamp_pdf()
#
# Returns:
#   seurat_obj — merged Seurat object with DecontX metadata and sample covariates
################################################################################

run_decontx <- function() {

  message("\n── DecontX: ambient RNA removal ─────────────────────────────────────────────")

  samples <- SAMPLE_METADATA$sample_id

  # ── Read counts and run DecontX ──────────────────────────────────────────────
  # Samples are independent — run in parallel across N_WORKERS

  message("  Running DecontX on ", length(samples), " sample(s) with ",
          N_WORKERS, " worker(s)...")

  sce_results <- future.apply::future_lapply(
    samples,
    function(s) {
      paths   <- cellranger_paths(s)
      sce     <- SingleCellExperiment(list(counts = Read10X(paths$filtered)))
      sce_raw <- SingleCellExperiment(list(counts = Read10X(paths$raw)))
      sce     <- decontX(sce, background = sce_raw)
      sce
    },
    future.seed = TRUE
  )
  sce_list <- setNames(sce_results, samples)

  # ── Save per-sample SCE objects and plot contamination UMAPs ─────────────────

  for (s in samples) {
    sce_obj <- sce_list[[s]]
    save(sce_obj, file = stamp("DecontX_SCE", s, ext = "RData"))
  }

  pdf(stamp_pdf("DecontX_contamination_UMAPs"), height = 5, width = 7)
  for (s in samples) {
    print(plotDecontXContamination(sce_list[[s]]) + ggtitle(s))
  }
  dev.off()
  message("  DecontX contamination UMAPs saved.")

  # ── Build per-sample Seurat objects from DecontX-corrected counts ─────────────
  # round() can push very low-count cells to all-zero UMI — drop them immediately.

  seurat_results <- future.apply::future_lapply(
    samples,
    function(s) {
      obj      <- CreateSeuratObject(
        counts    = round(decontXcounts(sce_list[[s]])),
        meta.data = as.data.frame(colData(sce_list[[s]]))
      )
      n_before <- ncol(obj)
      obj      <- obj[, Matrix::colSums(obj@assays$RNA@counts) > 0]
      n_drop   <- n_before - ncol(obj)
      if (n_drop > 0)
        message("  [", s, "] Dropped ", n_drop,
                " cell(s) with zero UMI after DecontX rounding")
      obj
    },
    future.seed = TRUE
  )
  seurat_list <- setNames(seurat_results, samples)
  rm(sce_list)

  # ── Merge all samples ────────────────────────────────────────────────────────

  if (length(seurat_list) == 1) {
    seurat_obj <- seurat_list[[1]]
  } else {
    seurat_obj <- merge(
      seurat_list[[1]],
      y            = seurat_list[-1],
      add.cell.ids = samples,
      project      = PROJECT_NAME
    )
  }

  rm(sce_list, seurat_list)
  gc()

  # ── Attach sample-level metadata to cells ────────────────────────────────────
  # Matches sample_id prefix in cell barcodes added by merge()

  cell_names <- colnames(seurat_obj)
  for (col in setdiff(colnames(SAMPLE_METADATA), "fastq_sample_name")) {
    meta_vec <- rep(NA_character_, length(cell_names))
    for (i in seq_len(nrow(SAMPLE_METADATA))) {
      idx            <- grep(paste0("^", SAMPLE_METADATA$sample_id[i], "_"), cell_names)
      meta_vec[idx]  <- as.character(SAMPLE_METADATA[[col]][i])
    }
    seurat_obj <- AddMetaData(seurat_obj, metadata = meta_vec, col.name = col)
  }

  save(seurat_obj, file = stamp("post_DecontX_Seurat", ext = "RData"))
  message("Post-DecontX merged Seurat object saved.")

  invisible(seurat_obj)
}
