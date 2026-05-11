################################################################################
# config.R — Project configuration
#
# Shared by both prepare_cellranger.R and pipeline.R.
# Edit this file and samples.tsv for each new project.
# Do not edit prepare_cellranger.R or pipeline.R.
#
# Workflow:
#   1. Fill in this file and samples.tsv
#   2. Rscript prepare_cellranger.R config.R               → validate FASTQs + generate shell script
#   3. Submit shell script on HPC; wait for jobs to finish
#   4. Rscript pipeline.R config.R                     → extract QC summary + run preprocessing
################################################################################

# ── Shared ────────────────────────────────────────────────────────────────────

# Project name — used as a prefix in all output file names
PROJECT_NAME <- "my_project"

# Directory where all output files (RData, PDFs, txt) will be written
OUTPUT_DIR <- "."

# Path to sample metadata file (.tsv or .csv); see samples.tsv for format
#
# Required columns:
#   sample_id         — unique sample label; used as cell-barcode prefix
#   fastq_sample_name — the --sample argument passed to CellRanger
#
# Optional columns (add or rename freely):
#   Any additional columns are automatically attached as per-cell metadata
#   in the Seurat object AND used as covariates in SCTransform regression.
#   Column names become metadata slot names — no hardcoding needed elsewhere.
METADATA_FILE <- "samples.tsv"   # use "samples.csv" if comma-delimited

# ── CellRanger (prepare_cellranger.R) ─────────────────────────────────────────────

CELLRANGER_PATH    <- "/path/to/cellranger"
TRANSCRIPTOME_PATH <- "/path/to/cellranger/reference"
FASTQ_DIR          <- "/path/to/fastq/files"
CELLRANGER_OUTDIR  <- "/path/to/cellranger/output"   # where CellRanger writes results
LOCALMEM           <- 64    # GB of RAM to allocate per job
LOCALCORES         <- 12    # CPU cores per job

# ── Preprocessing (pipeline.R) ────────────────────────────────────────────────

# Mitochondrial gene detection
# MITO_PATTERN: regex to identify mitochondrial genes
#   "^mt-"  for mouse  (e.g. mt-Nd1, mt-Co1)
#   "^MT-"  for human  (e.g. MT-ND1, MT-CO1)
MITO_PATTERN <- "^mt-"

# MITO_COL: name of the metadata column Seurat will store the mito % in.
# Convention: "percent.mito" (mouse) or "percent.mt" (human) but any name works.
# Must match what you use in SCT_EXTRA_VARS below.
MITO_COL     <- "percent.mito"

# QC thresholds
QC_MIN_FEATURES <- 500    # minimum genes detected per cell
QC_MAX_FEATURES <- 1e5    # maximum genes detected per cell (removes likely multiplets)
QC_MAX_MITO_PCT <- 15     # maximum mitochondrial gene %
QC_MAX_DECONTX  <- 0.25   # maximum DecontX contamination score
QC_MIN_CELLS    <- 20     # minimum cells a gene must appear in to be retained

# Clustering resolution
# CLUSTER_RESOLUTION: set to a numeric value to use a fixed resolution, or
#                     set to NULL to automatically select using clustree.
#                     Typical range: 0.2 (fewer, broader clusters) to 2.0
#                     (many fine-grained clusters). Higher values = more clusters.
CLUSTER_RESOLUTION <- NULL       # e.g. 0.8, or NULL for clustree auto-selection

# CLUSTREE_RESOLUTIONS: the range of resolutions tested when CLUSTER_RESOLUTION
#                       is NULL. Clustree picks the resolution where the tree
#                       first becomes stable (no new clusters forming).
CLUSTREE_RESOLUTIONS <- seq(0.2, 2.0, by = 0.2)

# Dimensionality reduction
N_PCS <- 50   # number of PCs to compute; optimal cutoff selected automatically

# UMAP settings
# umap.method: "uwot" (R-native, default) or "umap-learn" (Python via reticulate)
# metric:      "cosine" (recommended for scRNA-seq) or "correlation", "euclidean"
UMAP_METHOD <- "uwot"
UMAP_METRIC <- "cosine"

# SCTransform regression variables
# Only technical/unwanted sources of variation should be regressed out.
# Do NOT include biological variables (age, treatment, genotype) here —
# regressing those out removes the signal you want to study.
# decontX_contamination is added automatically.
SCT_EXTRA_VARS <- c("nFeature_RNA", MITO_COL)

# BATCH_COLS: technical batch column(s) to regress out and evaluate.
# Set to character(0) if there are no batch variables (default).
# Must exactly match column names in METADATA_FILE.
BATCH_COLS <- character(0)           # e.g. c("batch") or c("batch", "run")

# Doublet detection method
# Controls which doublet caller(s) are used and how calls are combined.
#   "union"         — flag as doublet if called by DoubletFinder OR scds (default, most conservative)
#   "intersection"  — flag as doublet only if called by BOTH methods (fewer removals)
#   "DoubletFinder" — use DoubletFinder result only
#   "scds"          — use scds result only
DOUBLET_METHOD <- "union"

# Extra doublet rate buffer on top of the 10x-predicted rate to account
DOUBLET_EXTRA_RATE <- 0.025
