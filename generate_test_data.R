################################################################################
# generate_test_data.R
#
# Creates a self-contained test directory with:
#   - Synthetic CellRanger-style sparse matrix output (filtered + raw)
#   - Synthetic metrics_summary.csv per sample
#   - FASTQ stub files in three naming scenarios:
#       scenario_a: all files correctly named  → expect 0 issues
#       scenario_b: mixed naming issues        → expect flagging + rename
#       scenario_c: files missing entirely     → expect hard stop
#   - Three matching config + samples files for each scenario
#
# Usage:
#   Rscript generate_test_data.R
#   # Then test each scenario:
#   Rscript prepare_cellranger.R test/scenario_a/config.R
#   Rscript prepare_cellranger.R test/scenario_b/config.R
#   Rscript prepare_cellranger.R test/scenario_b/config.R --rename
#   Rscript prepare_cellranger.R test/scenario_c/config.R
#   Rscript pipeline.R       test/scenario_a/config.R
################################################################################

library(Matrix)

set.seed(42)

TEST_ROOT <- "test"

# ── Parameters ────────────────────────────────────────────────────────────────

SAMPLES <- c("SampleA", "SampleB")  # fastq_sample_name values

SAMPLE_IDS <- c("CTL_1", "CTL_2")   # sample_id values

N_GENES_FILTERED <- 200   # genes in filtered matrix
N_CELLS_FILTERED <- 100   # cells per sample in filtered matrix
N_GENES_RAW      <- 200   # same gene set for simplicity
N_CELLS_RAW      <- 500   # larger empty-droplet pool for raw matrix

# Gene names: mix of mito and regular (mouse-style)
GENE_NAMES <- c(
  paste0("mt-Gene", seq_len(20)),          # mitochondrial
  paste0("Gene", seq_len(N_GENES_FILTERED - 20))
)

################################################################################
# HELPER: write a 10x-style sparse matrix directory
################################################################################

write_10x_dir <- function(dir_path, count_mat, gene_names, barcode_names) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  # barcodes.tsv.gz
  gz <- gzcon(file(file.path(dir_path, "barcodes.tsv.gz"), "wb"))
  writeLines(barcode_names, gz)
  close(gz)

  # features.tsv.gz  (gene_id, gene_name, feature_type)
  feat_df <- data.frame(
    id   = paste0("ENSMUSG", formatC(seq_along(gene_names), width = 11, flag = "0")),
    name = gene_names,
    type = "Gene Expression"
  )
  gz <- gzcon(file(file.path(dir_path, "features.tsv.gz"), "wb"))
  writeLines(apply(feat_df, 1, paste, collapse = "\t"), gz)
  close(gz)

  # matrix.mtx.gz
  sp   <- as(count_mat, "dgCMatrix")
  tmp  <- tempfile(fileext = ".mtx")
  writeMM(sp, tmp)
  raw_lines <- readLines(tmp)
  # Inject MatrixMarket header comment expected by Read10X
  raw_lines[1] <- "%%MatrixMarket matrix coordinate integer general"
  gz <- gzcon(file(file.path(dir_path, "matrix.mtx.gz"), "wb"))
  writeLines(raw_lines, gz)
  close(gz)
}

################################################################################
# HELPER: write a realistic metrics_summary.csv
################################################################################

write_metrics_csv <- function(path, n_cells) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  metrics <- data.frame(
    `Estimated Number of Cells`            = n_cells,
    `Mean Reads per Cell`                  = sample(25000:50000, 1),
    `Median Genes per Cell`                = sample(1500:3000,   1),
    `Number of Reads`                      = n_cells * sample(25000:50000, 1),
    `Valid Barcodes`                       = paste0(round(runif(1, 95, 99), 1), "%"),
    `Sequencing Saturation`                = paste0(round(runif(1, 40, 80), 1), "%"),
    `Q30 Bases in Barcode`                 = paste0(round(runif(1, 93, 98), 1), "%"),
    `Q30 Bases in RNA Read`                = paste0(round(runif(1, 88, 95), 1), "%"),
    `Q30 Bases in UMI`                     = paste0(round(runif(1, 93, 98), 1), "%"),
    `Reads Mapped to Genome`               = paste0(round(runif(1, 88, 96), 1), "%"),
    `Reads Mapped Confidently to Genome`   = paste0(round(runif(1, 80, 92), 1), "%"),
    `Reads Mapped Confidently to Transcriptome` = paste0(round(runif(1, 55, 75), 1), "%"),
    `Median UMI Counts per Cell`           = sample(3000:8000, 1),
    check.names = FALSE
  )
  write.csv(metrics, file = path, row.names = FALSE)
}

################################################################################
# HELPER: generate synthetic count matrix
################################################################################

make_count_matrix <- function(n_genes, n_cells, sparsity = 0.85) {
  mat <- matrix(
    rbinom(n_genes * n_cells, size = 20, prob = 0.3) *
      rbinom(n_genes * n_cells, size = 1, prob = 1 - sparsity),
    nrow = n_genes, ncol = n_cells
  )
  # Ensure mito genes (first 20) have somewhat elevated counts in some cells
  mat[1:20, sample(n_cells, n_cells %/% 5)] <- mat[1:20, sample(n_cells, n_cells %/% 5)] + 
    matrix(rpois(20 * (n_cells %/% 5), lambda = 5), nrow = 20)
  mat
}

################################################################################
# GENERATE CELLRANGER OUTPUT (shared across all scenarios)
################################################################################

cr_out_dir <- file.path(TEST_ROOT, "cellranger_output")
message("Generating synthetic CellRanger output in: ", cr_out_dir)

for (i in seq_along(SAMPLE_IDS)) {
  sid <- SAMPLE_IDS[i]
  message("  Sample: ", sid)

  # ── Filtered matrix ──────────────────────────────────────────────────────────
  barcodes_filt <- paste0("CELL", formatC(seq_len(N_CELLS_FILTERED), width = 6, flag = "0"), "-1")
  counts_filt   <- make_count_matrix(N_GENES_FILTERED, N_CELLS_FILTERED, sparsity = 0.75)
  rownames(counts_filt) <- GENE_NAMES
  colnames(counts_filt) <- barcodes_filt

  filt_dir <- file.path(cr_out_dir, sid, "outs", "filtered_feature_bc_matrix")
  write_10x_dir(filt_dir, counts_filt, GENE_NAMES, barcodes_filt)

  # ── Raw matrix (filtered cells + empty droplets) ─────────────────────────────
  barcodes_raw <- c(
    barcodes_filt,
    paste0("EMPTY", formatC(seq_len(N_CELLS_RAW - N_CELLS_FILTERED), width = 6, flag = "0"), "-1")
  )
  counts_raw_extra <- make_count_matrix(N_GENES_FILTERED, N_CELLS_RAW - N_CELLS_FILTERED,
                                         sparsity = 0.98)  # very sparse empty droplets
  counts_raw <- cbind(counts_filt, counts_raw_extra)
  colnames(counts_raw) <- barcodes_raw

  raw_dir <- file.path(cr_out_dir, sid, "outs", "raw_feature_bc_matrix")
  write_10x_dir(raw_dir, counts_raw, GENE_NAMES, barcodes_raw)

  # ── metrics_summary.csv ───────────────────────────────────────────────────────
  write_metrics_csv(
    file.path(cr_out_dir, sid, "outs", "metrics_summary.csv"),
    N_CELLS_FILTERED
  )
}

message("CellRanger output written.\n")

################################################################################
# SCENARIO A: All FASTQs correctly named — expect 0 issues
################################################################################

message("── Scenario A: correctly named FASTQs ────────────────────────────────────────")

fastq_dir_a <- file.path(TEST_ROOT, "scenario_a", "fastqs")
dir.create(fastq_dir_a, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(SAMPLES)) {
  sname <- SAMPLES[i]
  for (lane in c("L001", "L002")) {
    for (read in c("R1", "R2")) {
      fname <- sprintf("%s_S%d_%s_%s_001.fastq.gz", sname, i, lane, read)
      file.create(file.path(fastq_dir_a, fname))
    }
  }
}

message("  Files: ", paste(list.files(fastq_dir_a), collapse = ", "))

# config
writeLines(sprintf('
PROJECT_NAME    <- "test_scenario_a"
OUTPUT_DIR      <- "test/scenario_a/output"
METADATA_FILE   <- "test/scenario_a/samples.tsv"

CELLRANGER_PATH    <- "/path/to/cellranger"
TRANSCRIPTOME_PATH <- "/path/to/reference"
FASTQ_DIR          <- "%s"
CELLRANGER_OUTDIR  <- "%s"
LOCALMEM           <- 8
LOCALCORES         <- 2

MITO_PATTERN       <- "^mt-"
MITO_COL           <- "percent.mito"
QC_MIN_FEATURES    <- 50
QC_MAX_FEATURES    <- 1e5
QC_MAX_MITO_PCT    <- 50
QC_MAX_DECONTX     <- 0.9
QC_MIN_CELLS       <- 2
N_PCS              <- 10
SCT_EXTRA_VARS     <- c("nFeature_RNA", MITO_COL)
BATCH_COLS         <- c("batch")
DOUBLET_METHOD     <- "union"
DOUBLET_EXTRA_RATE <- 0.025
', fastq_dir_a, cr_out_dir),
  file.path(TEST_ROOT, "scenario_a", "config.R")
)

# samples.tsv
write.table(
  data.frame(
    sample_id         = SAMPLE_IDS,
    fastq_sample_name = SAMPLES,
    condition         = c("CTL", "CTL"),
    batch             = c("batch1", "batch2")
  ),
  file      = file.path(TEST_ROOT, "scenario_a", "samples.tsv"),
  sep       = "\t", row.names = FALSE, quote = FALSE
)

dir.create(file.path(TEST_ROOT, "scenario_a", "output"), recursive = TRUE, showWarnings = FALSE)
message("  Scenario A config written.\n")

################################################################################
# SCENARIO B: Mixed naming issues — expect report + successful --rename
################################################################################

message("── Scenario B: mixed naming issues ──────────────────────────────────────────")

fastq_dir_b <- file.path(TEST_ROOT, "scenario_b", "fastqs")
dir.create(fastq_dir_b, recursive = TRUE, showWarnings = FALSE)

# SampleA: correct
file.create(file.path(fastq_dir_b, "SampleA_S1_L001_R1_001.fastq.gz"))
file.create(file.path(fastq_dir_b, "SampleA_S1_L001_R2_001.fastq.gz"))

# SampleB: various issues
file.create(file.path(fastq_dir_b, "SampleB_lane1_R1.fastq.gz"))      # missing S-index, wrong lane format
file.create(file.path(fastq_dir_b, "SampleB_S1_L001_R2.fastq.gz"))    # missing _001 suffix
file.create(file.path(fastq_dir_b, "SampleB_S2_L002_read1_001.fastq.gz")) # "read1" instead of "R1"

message("  Files: ", paste(list.files(fastq_dir_b), collapse = ", "))

writeLines(sprintf('
PROJECT_NAME    <- "test_scenario_b"
OUTPUT_DIR      <- "test/scenario_b/output"
METADATA_FILE   <- "test/scenario_b/samples.tsv"

CELLRANGER_PATH    <- "/path/to/cellranger"
TRANSCRIPTOME_PATH <- "/path/to/reference"
FASTQ_DIR          <- "%s"
CELLRANGER_OUTDIR  <- "%s"
LOCALMEM           <- 8
LOCALCORES         <- 2

MITO_PATTERN       <- "^mt-"
MITO_COL           <- "percent.mito"
QC_MIN_FEATURES    <- 50
QC_MAX_FEATURES    <- 1e5
QC_MAX_MITO_PCT    <- 50
QC_MAX_DECONTX     <- 0.9
QC_MIN_CELLS       <- 2
N_PCS              <- 10
SCT_EXTRA_VARS     <- c("nFeature_RNA", MITO_COL)
BATCH_COLS         <- c("batch")
DOUBLET_METHOD     <- "union"
DOUBLET_EXTRA_RATE <- 0.025
', fastq_dir_b, cr_out_dir),
  file.path(TEST_ROOT, "scenario_b", "config.R")
)

write.table(
  data.frame(
    sample_id         = SAMPLE_IDS,
    fastq_sample_name = SAMPLES,
    condition         = c("CTL", "CTL"),
    batch             = c("batch1", "batch2")
  ),
  file      = file.path(TEST_ROOT, "scenario_b", "samples.tsv"),
  sep       = "\t", row.names = FALSE, quote = FALSE
)

dir.create(file.path(TEST_ROOT, "scenario_b", "output"), recursive = TRUE, showWarnings = FALSE)
message("  Scenario B config written.\n")

################################################################################
# SCENARIO C: FASTQ files missing entirely — expect hard stop
################################################################################

message("── Scenario C: missing FASTQs ────────────────────────────────────────────────")

fastq_dir_c <- file.path(TEST_ROOT, "scenario_c", "fastqs")
dir.create(fastq_dir_c, recursive = TRUE, showWarnings = FALSE)
# Intentionally empty — no files created

writeLines(sprintf('
PROJECT_NAME    <- "test_scenario_c"
OUTPUT_DIR      <- "test/scenario_c/output"
METADATA_FILE   <- "test/scenario_c/samples.tsv"

CELLRANGER_PATH    <- "/path/to/cellranger"
TRANSCRIPTOME_PATH <- "/path/to/reference"
FASTQ_DIR          <- "%s"
CELLRANGER_OUTDIR  <- "%s"
LOCALMEM           <- 8
LOCALCORES         <- 2

MITO_PATTERN       <- "^mt-"
MITO_COL           <- "percent.mito"
QC_MIN_FEATURES    <- 50
QC_MAX_FEATURES    <- 1e5
QC_MAX_MITO_PCT    <- 50
QC_MAX_DECONTX     <- 0.9
QC_MIN_CELLS       <- 2
N_PCS              <- 10
SCT_EXTRA_VARS     <- c("nFeature_RNA", MITO_COL)
BATCH_COLS         <- c("batch")
DOUBLET_METHOD     <- "union"
DOUBLET_EXTRA_RATE <- 0.025
', fastq_dir_c, cr_out_dir),
  file.path(TEST_ROOT, "scenario_c", "config.R")
)

write.table(
  data.frame(
    sample_id         = SAMPLE_IDS,
    fastq_sample_name = SAMPLES,
    condition         = c("CTL", "CTL"),
    batch             = c("batch1", "batch2")
  ),
  file      = file.path(TEST_ROOT, "scenario_c", "samples.tsv"),
  sep       = "\t", row.names = FALSE, quote = FALSE
)

dir.create(file.path(TEST_ROOT, "scenario_c", "output"), recursive = TRUE, showWarnings = FALSE)
message("  Scenario C config written.\n")

################################################################################
# SUMMARY
################################################################################

message("════════════════════════════════════════════════════════════════════════════════")
message("Test data generated under: ", normalizePath(TEST_ROOT))
message("")
message("Test prepare_cellranger.R:")
message("  Scenario A (valid names)  → Rscript prepare_cellranger.R test/scenario_a/config.R")
message("             expected: 0 issues, shell script generated")
message("")
message("  Scenario B (bad names)    → Rscript prepare_cellranger.R test/scenario_b/config.R")
message("             expected: issues flagged, report written")
message("           + --rename       → Rscript prepare_cellranger.R test/scenario_b/config.R --rename")
message("             expected: files renamed, rename log written")
message("")
message("  Scenario C (no FASTQs)    → Rscript prepare_cellranger.R test/scenario_c/config.R")
message("             expected: 'No FASTQ files found' for all samples, warning printed")
message("")
message("Test pipeline.R (uses shared CellRanger output, relaxed QC thresholds):")
message("  Rscript pipeline.R test/scenario_a/config.R")
message("════════════════════════════════════════════════════════════════════════════════")
