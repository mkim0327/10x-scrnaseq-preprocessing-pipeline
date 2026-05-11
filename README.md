# scRNA-seq Preprocessing Pipeline

A generalizable R pipeline for preprocessing 10x Genomics single-cell RNA-seq data, covering FASTQ validation, CellRanger job generation, ambient RNA removal, doublet detection, and Seurat QC filtering.

> **Note:** This pipeline was developed with assistance from [Claude](https://claude.ai) (Anthropic).

---

## Overview

The pipeline is split into two independent scripts that reflect the natural boundary between HPC compute (CellRanger) and local/interactive R analysis (preprocessing):

```
prepare_cellranger.R   →   [HPC: run_cellranger.sh]   →   pipeline.R
```

| Script | What it does |
|--------|-------------|
| `prepare_cellranger.R` | Validates FASTQ naming, flags/renames non-conforming files, generates a CellRanger shell script for HPC submission |
| `pipeline.R` | Reads CellRanger output; runs DecontX, doublet detection, QC visualization, and Seurat filtering |

Both scripts are fully driven by `config.R` and a sample metadata file. 

---

## Repository Structure

```
├── prepare_cellranger.R     # FASTQ validation + CellRanger script generator
├── pipeline.R               # Full R preprocessing pipeline
├── config.R                 # Project configuration template (edit per project)
├── samples_cellranger.tsv   # Sample sheet template
│
├── config_cellranger.R      # Example filled-in config for prepare_cellranger.R
├── samples_VCD_CTL.tsv      # Example sample sheet (Benayoun lab VCD CTL dataset)
│
└── tests/
    ├── generate_test_data.R     # Generates synthetic CellRanger output + FASTQ scenarios
    ├── config_test_cellranger.R # Test config for prepare_cellranger.R
    └── config_test_pipeline.R   # Test config for pipeline.R (relaxed QC thresholds)
```

---

## Prerequisites

### CellRanger
- [CellRanger](https://www.10xgenomics.com/support/software/cell-ranger) ≥ 6.0
- A matching genome reference (e.g. mouse: `refdata-gex-mm10-2020-A`, human: `refdata-gex-GRCh38-2020-A`)

### R packages

| Package | Source | Purpose |
|---------|--------|---------|
| `Seurat` | CRAN | Single-cell analysis framework |
| `sctransform` | CRAN | Variance-stabilizing normalization |
| `celda` | Bioconductor | DecontX ambient RNA removal |
| `singleCellTK` | Bioconductor | SCE utilities |
| `SingleCellExperiment` | Bioconductor | Data container |
| `DoubletFinder` | GitHub | Doublet detection |
| `scds` | Bioconductor | Hybrid doublet scoring |
| `scater` | Bioconductor | QC utilities |
| `ggplot2` | CRAN | Visualization |
| `reshape2` | CRAN | Data reshaping for plots |
| `dplyr` | CRAN | Data manipulation |
| `clustree` | CRAN | Clustering resolution selection |
| `bitops` | CRAN | Bitwise operations |

Install from Bioconductor:
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("celda", "singleCellTK", "SingleCellExperiment", "scds", "scater"))
```

Install DoubletFinder from GitHub:
```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

Install remaining packages from CRAN:
```r
install.packages(c("Seurat", "sctransform", "ggplot2", "reshape2",
                   "dplyr", "clustree", "scales", "bitops", "reshape2"))
```

---

## Quick Start

### 1. Set up your project files

Copy and fill in the templates:

```bash
cp config.R             my_project_config.R
cp samples_cellranger.tsv  my_project_samples.tsv
```

Edit `my_project_config.R` — all fields marked `<FILL IN>` must be set. Edit `my_project_samples.tsv` — one row per sample (see [Sample Sheet](#sample-sheet) below).

### 2. Validate FASTQs and generate CellRanger script

```bash
Rscript prepare_cellranger.R my_project_config.R
```

If naming issues are found, review the report and re-run with `--rename` to fix automatically:

```bash
Rscript prepare_cellranger.R my_project_config.R --rename
```

This writes a shell script `YYYY-MM-DD_<PROJECT>_run_cellranger.sh` to `OUTPUT_DIR`.

### 3. Submit CellRanger jobs on HPC

```bash
bash YYYY-MM-DD_my_project_run_cellranger.sh
```

Or adapt to your scheduler (SLURM, SGE, etc.) as needed.

### 4. Run the preprocessing pipeline

Once all CellRanger jobs are complete:

```bash
Rscript pipeline.R my_project_config.R
```

---

## Pipeline Steps

### `prepare_cellranger.R`

**Step 1 — FASTQ validation**
Scans `FASTQ_DIR` for all `.fastq.gz` files and checks each against the CellRanger required naming convention:
```
[SampleName]_S[Index]_L00[Lane]_[ReadType]_001.fastq.gz
```
Issues are written to a TSV report. With `--rename`, parseable filenames are renamed in place and a rename log is saved.

**Step 2 — Shell script generation**
Writes one `cellranger count` command per sample, ready to submit to an HPC cluster.

---

### `pipeline.R`

**Step 1 — CellRanger experiment summary**
Reads `metrics_summary.csv` from each sample's CellRanger output and merges into a single summary table. Produces:
- Bar charts per metric per sample, coloured by each metadata covariate
- Z-scored heatmap across all metrics and samples (outlier detection)
- Grouped chart of all sequencing quality percentage metrics

**Step 2 — Ambient RNA removal (DecontX)**
Runs `decontX()` on each sample using the raw (unfiltered) matrix as background. Saves per-sample SCE objects and contamination UMAP plots.

**Step 2b — QC visualizations**
Runs on the full unfiltered merged object before any cell removal. Produces:
- Violin + boxplot for each QC metric × each metadata grouping variable
- DecontX contamination density curves per group (threshold line overlaid)
- Pairwise scatter plots coloured by DecontX score with QC threshold lines
- Cell count bar charts per sample and per covariate
- Pearson correlation heatmap between QC metrics

**Step 3 — Doublet detection**
Doublets are identified *before* QC filtering to preserve the cell context that doublet callers rely on. Uses a union call of two methods:
- **DoubletFinder** — pN/pK sweep with per-sample doublet rate estimated from the 10x Genomics cell recovery table
- **scds (cxds_bcds_hybrid)** — co-expression and count-based hybrid scoring

A cell is called a doublet if flagged by *either* method.

**Step 4 — Seurat QC filtering**
Filters cells by:

| Parameter | Config variable | Default |
|-----------|----------------|---------|
| Min genes per cell | `QC_MIN_FEATURES` | 500 |
| Max genes per cell | `QC_MAX_FEATURES` | 100,000 |
| Max mitochondrial % | `QC_MAX_MITO_PCT` | 15 |
| Max DecontX score | `QC_MAX_DECONTX` | 0.25 |

Followed by SCTransform normalization, PCA, UMAP, and Louvain clustering. PC cutoff is selected automatically using cumulative variance (>90%) and consecutive variance drop (<0.1%) criteria.

---

## Configuration

All project-specific settings live in `config.R`. The scripts never need to be edited.

### Key parameters

```r
PROJECT_NAME  <- "my_project"          # prefix for all output files
OUTPUT_DIR    <- "/path/to/output"     # data files written here
METADATA_FILE <- "samples.tsv"        # path to sample sheet

# Species
MITO_PATTERN  <- "^mt-"               # "^mt-" mouse / "^MT-" human
MITO_COL      <- "percent.mito"       # "percent.mito" mouse / "percent.mt" human
```

### Sample sheet

The metadata file requires exactly two columns; all additional columns are automatically used as per-cell metadata in Seurat *and* as SCTransform regression covariates — no further configuration needed.

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | ✓ | Unique sample label; used as CellRanger `--id` and cell barcode prefix |
| `fastq_sample_name` | ✓ | Sample name in FASTQ filenames; passed to CellRanger `--sample` |
| `age`, `treatment`, `batch`, ... | optional | Any experimental grouping variables |

Example:
```
sample_id         fastq_sample_name  age         batch
Young_female_1    MHK305             4m          batch1
Young_female_2    MHK306             4m          batch2
```

### Output directory layout

```
OUTPUT_DIR/
├── plots/                          # all PDF figures
│   ├── YYYY-MM-DD_<PROJECT>_CellRanger_summary_by_sample_id.pdf
│   ├── YYYY-MM-DD_<PROJECT>_CellRanger_summary_heatmap.pdf
│   ├── YYYY-MM-DD_<PROJECT>_CellRanger_sequencing_quality.pdf
│   ├── YYYY-MM-DD_<PROJECT>_DecontX_contamination_UMAPs.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_violin_<group>.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_decontX_density_<group>.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_scatter_metrics.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_cell_counts.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_metric_correlation_heatmap.pdf
│   ├── YYYY-MM-DD_<PROJECT>_DoubletFinder_UMAPs.pdf
│   ├── YYYY-MM-DD_<PROJECT>_scds_UMAPs.pdf
│   ├── YYYY-MM-DD_<PROJECT>_DoubletCall_UMAP.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_violinPlots_PRE_FILTER.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_scatter_PRE_FILTER.pdf
│   ├── YYYY-MM-DD_<PROJECT>_QC_violinPlots_POST_FILTER.pdf
│   ├── YYYY-MM-DD_<PROJECT>_ElbowPlot_threshold.pdf
│   └── YYYY-MM-DD_<PROJECT>_UMAP_SeuratClusters_res2.0.pdf
│
├── YYYY-MM-DD_<PROJECT>_CellRanger_QC_summary.txt
├── YYYY-MM-DD_<PROJECT>_DecontX_SCE_<sample>.RData
├── YYYY-MM-DD_<PROJECT>_post_DecontX_Seurat.RData
├── YYYY-MM-DD_<PROJECT>_Seurat_AnnotatedDoublets.RData
├── YYYY-MM-DD_<PROJECT>_Seurat_post_QC.RData
├── YYYY-MM-DD_<PROJECT>_Seurat_final.RData
└── YYYY-MM-DD_<PROJECT>_session_info.txt
```

---

## Testing

Generate synthetic CellRanger output and test FASTQ scenarios:

```bash
Rscript tests/generate_test_data.R
```

This creates three scenarios under `test/`:

| Scenario | FASTQ files | Expected result |
|----------|-------------|-----------------|
| `scenario_a` | All correctly named | 0 issues, shell script generated |
| `scenario_b` | Mixed naming problems | Issues flagged; `--rename` fixes them |
| `scenario_c` | No files found | Warning for all samples |

Test each:
```bash
Rscript prepare_cellranger.R test/scenario_a/config.R
Rscript prepare_cellranger.R test/scenario_b/config.R
Rscript prepare_cellranger.R test/scenario_b/config.R --rename
Rscript prepare_cellranger.R test/scenario_c/config.R

# Full pipeline test (uses synthetic CellRanger output)
Rscript pipeline.R test/scenario_a/config.R
```

---

## Notes

- **SCTransform only regresses out technical variables**, not all metadata columns. Biological variables (`age`, `genotype`, `treatment`, etc.) are attached as per-cell metadata and excluded from regression. Only the variables listed in `BATCH_COLS` (plus `nFeature_RNA`, `MITO_COL`, and `decontX_contamination`) are regressed out. To add a new technical covariate, add its column name to `BATCH_COLS` in `config.R`.
- **Species switching** — to use with human data, set `MITO_PATTERN <- "^MT-"` and `MITO_COL <- "percent.mt"` in `config.R`. Nothing else needs to change.

---

## Acknowledgements

Pipeline code generated with assistance from [Claude](https://claude.ai) (claude-sonnet-4-6, Anthropic).

### References

- **CellRanger**: 10x Genomics. https://www.10xgenomics.com/support/software/cell-ranger
- **Seurat**: Hao et al. (2021). *Integrated analysis of multimodal single-cell data.* Cell. https://doi.org/10.1016/j.cell.2021.04.048
- **DecontX**: Yang et al. (2020). *Decontamination of ambient RNA in single-cell RNA-seq with DecontX.* Genome Biology. https://doi.org/10.1186/s13059-020-1950-6
- **DoubletFinder**: McGinnis et al. (2019). *DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors.* Cell Systems. https://doi.org/10.1016/j.cels.2019.03.003
- **scds**: Bais & Kostka (2020). *scds: Computational annotation of doublets in single-cell RNA sequencing data.* Bioinformatics. https://doi.org/10.1093/bioinformatics/btz698
- **SCTransform**: Hafemeister & Satija (2019). *Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.* Genome Biology. https://doi.org/10.1186/s13059-019-1874-1
