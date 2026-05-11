################################################################################
# modules/utils.R — Shared helper functions
################################################################################

#' Build CellRanger output paths for a sample
cellranger_paths <- function(sid) {
  base <- file.path(CELLRANGER_OUTDIR, sid, "outs")
  list(
    filtered = file.path(base, "filtered_feature_bc_matrix"),
    raw      = file.path(base, "raw_feature_bc_matrix"),
    metrics  = file.path(base, "metrics_summary.csv")
  )
}

#' Auto-select number of PCs:
#'   co1 = first PC where cumulative variance > 90% AND per-PC variance < 5%
#'   co2 = last PC where consecutive drop in variance > 0.1%
#'   returns min(co1, co2); falls back to npcs/2 if either criterion cannot be met
select_pcs <- function(seurat_obj, reduction = "pca") {
  pct  <- seurat_obj[[reduction]]@stdev / sum(seurat_obj[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  npcs <- length(pct)

  co1_idx <- which(cumu > 90 & pct < 5)
  co1     <- if (length(co1_idx) > 0) co1_idx[1] else npcs

  co2_idx <- sort(which((pct[1:(npcs - 1)] - pct[2:npcs]) > 0.1), decreasing = TRUE)
  co2     <- if (length(co2_idx) > 0) co2_idx[1] + 1 else npcs

  pcs <- min(co1, co2, npcs)
  message("  PC selection: co1=", co1, "  co2=", co2, "  using ", pcs, " PCs")
  pcs
}

#' Build a dated, project-stamped path for data output files (RData, txt, etc.)
stamp <- function(..., ext) {
  fname <- paste(c(as.character(Sys.Date()), PROJECT_NAME, ...), collapse = "_")
  file.path(OUTPUT_DIR, paste0(fname, ".", ext))
}

#' Build a dated, project-stamped path for PDF plot files
stamp_pdf <- function(...) {
  fname <- paste(c(as.character(Sys.Date()), PROJECT_NAME, ...), collapse = "_")
  file.path(PLOT_DIR, paste0(fname, ".pdf"))
}

#' Remove cells with zero total UMI — prevents log(0) = -Inf in SCTransform
drop_zero_umi <- function(obj, label = "") {
  keep   <- Matrix::colSums(obj@assays$RNA@counts) > 0
  n_drop <- sum(!keep)
  if (n_drop > 0)
    message(if (nchar(label) > 0) paste0("[", label, "] ") else "",
            "Dropping ", n_drop, " cell(s) with zero UMI before SCTransform")
  obj[, keep]
}

#' Filter SCTransform regression variables to those with >= 2 unique non-NA
#' values in the given object. Variables with only one level cannot be used
#' as contrasts in model.matrix and will cause SCTransform to error.
#' Returns NULL (not character(0)) when no valid variables remain, since
#' SCTransform expects NULL rather than an empty vector for no regression.
filter_regression_vars <- function(obj, vars) {
  if (length(vars) == 0) return(NULL)
  md    <- obj@meta.data
  valid <- vapply(vars, function(v) {
    if (!v %in% colnames(md)) return(FALSE)
    length(unique(md[[v]][!is.na(md[[v]])])) >= 2
  }, logical(1))
  dropped <- vars[!valid]
  if (length(dropped) > 0)
    message("  Dropping regression variable(s) with < 2 levels in this object: ",
            paste(dropped, collapse = ", "))
  result <- vars[valid]
  if (length(result) == 0) NULL else result
}

#' Linear model of 10x doublet rate vs recovered cell count
make_doublet_rate_model <- function() {
  tbl <- data.frame(
    cell_number = c(3000, 4000,  5000,  6000,  7000,  8000,
                    9000, 10000, 12000, 14000, 16000, 18000, 20000),
    dblt_rate   = c(2.4,  3.2,   4.0,   4.8,   5.6,   6.4,
                    7.2,  8.0,   9.6,  11.2,  12.8,  14.4,  16.0)
  )
  lm(dblt_rate ~ cell_number, data = tbl)
}
