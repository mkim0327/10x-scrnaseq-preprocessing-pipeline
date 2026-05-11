################################################################################
# config_test_pipeline.R
# Test config for pipeline.R
# Based on Benayoun lab VCD dataset (CTL samples only)
################################################################################

# ── Shared ────────────────────────────────────────────────────────────────────

PROJECT_NAME  <- "10x_ovary_Benayoun_lab_VCD"
OUTPUT_DIR    <- "/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Benayoun_lab_VCD_model/1_Preprocessing"
METADATA_FILE <- "samples_VCD_CTL.tsv"

# ── CellRanger output location (pipeline.R reads from here) ───────────────────

CELLRANGER_OUTDIR <- "/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD"

# ── Preprocessing ─────────────────────────────────────────────────────────────

MITO_PATTERN       <- "^mt-"
MITO_COL           <- "percent.mito"

QC_MIN_FEATURES    <- 500
QC_MAX_FEATURES    <- 1e5
QC_MAX_MITO_PCT    <- 15
QC_MAX_DECONTX     <- 0.25
QC_MIN_CELLS       <- 20

N_PCS              <- 50

SCT_EXTRA_VARS     <- c("nFeature_RNA", MITO_COL)
BATCH_COLS         <- c("batch")       # set to character(0) if no batches

CLUSTER_RESOLUTION   <- 0.5          # use fixed resolution for faster testing
CLUSTREE_RESOLUTIONS <- seq(0.2, 1.0, by = 0.2)
DOUBLET_METHOD       <- "union"
DOUBLET_EXTRA_RATE   <- 0.025
