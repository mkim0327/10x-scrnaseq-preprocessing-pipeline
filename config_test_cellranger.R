################################################################################
# config_test_cellranger.R
# Test config for prepare_cellranger.R
# Based on Benayoun lab VCD dataset (CTL samples only)
################################################################################

# ── Shared ────────────────────────────────────────────────────────────────────

PROJECT_NAME  <- "10x_ovary_Benayoun_lab_VCD"
OUTPUT_DIR    <- "/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Benayoun_lab_VCD_model/1_Preprocessing"
METADATA_FILE <- "samples_VCD_CTL.tsv"

# ── CellRanger ────────────────────────────────────────────────────────────────

CELLRANGER_PATH    <- "/home/minhooki/Softwares/cellranger-7.1.0/cellranger"
TRANSCRIPTOME_PATH <- "/home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A"
FASTQ_DIR          <- "/mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD"
CELLRANGER_OUTDIR  <- "/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD"
LOCALMEM           <- 64
LOCALCORES         <- 12
