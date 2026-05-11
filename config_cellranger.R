################################################################################
# config_cellranger.R — CellRanger configuration template
#
# Usage:
#   Rscript prepare_cellranger.R config_cellranger.R
#   Rscript prepare_cellranger.R config_cellranger.R --rename
#
# Fill in every field marked <FILL IN> before running.
################################################################################

# ── Shared ────────────────────────────────────────────────────────────────────

# Short project label used as a prefix in all output file names
PROJECT_NAME <- "<FILL IN>"                  # e.g. "10x_ovary_MyProject"

# Directory where output files (shell script, FASTQ report) will be written
OUTPUT_DIR   <- "<FILL IN>"                  # e.g. "/path/to/project/cellranger"

# Path to the sample metadata file (see samples_cellranger.tsv)
METADATA_FILE <- "samples_cellranger.tsv"

# ── CellRanger ────────────────────────────────────────────────────────────────

# Full path to the CellRanger executable
CELLRANGER_PATH    <- "<FILL IN>"            # e.g. "/home/user/software/cellranger-7.1.0/cellranger"

# Full path to the CellRanger genome reference directory
TRANSCRIPTOME_PATH <- "<FILL IN>"            # e.g. "/home/user/references/refdata-gex-mm10-2020-A"
                                             # Mouse: refdata-gex-mm10-2020-A
                                             # Human: refdata-gex-GRCh38-2020-A

# Directory containing the raw FASTQ files
FASTQ_DIR          <- "<FILL IN>"            # e.g. "/data/rawdata/MyProject"

# Directory where CellRanger will write its output (one subfolder per sample_id)
CELLRANGER_OUTDIR  <- "<FILL IN>"            # e.g. "/data/cellranger/MyProject"

# Computational resources
LOCALMEM   <- 64    # RAM in GB  — increase for larger datasets
LOCALCORES <- 12    # CPU cores  — should not exceed available cores on your HPC node
