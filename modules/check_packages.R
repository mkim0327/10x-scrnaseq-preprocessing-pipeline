################################################################################
# modules/check_packages.R — Package installation check
#
# Checks whether all required packages are installed. Missing packages are
# installed automatically from the correct source (CRAN, Bioconductor, GitHub).
# Already-installed packages are loaded without reinstalling.
#
# Call this before loading any libraries in pipeline.R or prepare_cellranger.R.
#
# Sources:
#   CRAN        — install.packages()
#   Bioconductor — BiocManager::install()
#   GitHub      — remotes::install_github()
################################################################################

check_and_install_packages <- function() {

  message("\n── Checking required packages ───────────────────────────────────────────────")

  # ── Package registry ─────────────────────────────────────────────────────────
  # Each entry: list(source, repo)
  #   source  — "CRAN", "Bioconductor", or "GitHub"
  #   repo    — for GitHub: "owner/repo"; for CRAN/Bioc: same as package name

  packages <- list(
    # CRAN
    Seurat            = list(source = "CRAN",         repo = "Seurat"),
    sctransform       = list(source = "CRAN",         repo = "sctransform"),
    clustree          = list(source = "CRAN",         repo = "clustree"),
    scales            = list(source = "CRAN",         repo = "scales"),
    dplyr             = list(source = "CRAN",         repo = "dplyr"),
    ggplot2           = list(source = "CRAN",         repo = "ggplot2"),
    reshape2          = list(source = "CRAN",         repo = "reshape2"),
    Matrix            = list(source = "CRAN",         repo = "Matrix"),
    remotes           = list(source = "CRAN",         repo = "remotes"),
    BiocManager       = list(source = "CRAN",         repo = "BiocManager"),
    bitops            = list(source = "CRAN",         repo = "bitops"),
    cluster           = list(source = "CRAN",         repo = "cluster"),
    future            = list(source = "CRAN",         repo = "future"),
    future.apply      = list(source = "CRAN",         repo = "future.apply"),
    furrr             = list(source = "CRAN",         repo = "furrr"),

    # Bioconductor
    celda             = list(source = "Bioconductor", repo = "celda"),
    singleCellTK      = list(source = "Bioconductor", repo = "singleCellTK"),
    scds              = list(source = "Bioconductor", repo = "scds"),
    scater            = list(source = "Bioconductor", repo = "scater"),
    SingleCellExperiment = list(source = "Bioconductor", repo = "SingleCellExperiment"),

    # GitHub
    DoubletFinder     = list(source = "GitHub",       repo = "chris-mcginnis-ucsf/DoubletFinder")
  )

  # ── Identify missing packages ─────────────────────────────────────────────────

  installed  <- rownames(installed.packages())
  pkg_names  <- names(packages)
  missing    <- pkg_names[!pkg_names %in% installed]

  if (length(missing) == 0) {
    message("All packages already installed.")
  } else {
    message("Missing packages: ", paste(missing, collapse = ", "))
  }

  # ── Install missing packages by source ────────────────────────────────────────

  if (length(missing) == 0) {
    message("Loading packages...")
    for (pkg in pkg_names)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message("All packages loaded successfully.\n")
    return(invisible(TRUE))
  }

  missing_cran <- missing[vapply(missing, function(p) packages[[p]]$source == "CRAN",         logical(1))]
  missing_bioc <- missing[vapply(missing, function(p) packages[[p]]$source == "Bioconductor", logical(1))]
  missing_gh   <- missing[vapply(missing, function(p) packages[[p]]$source == "GitHub",       logical(1))]

  # CRAN
  if (length(missing_cran) > 0) {
    message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  # Bioconductor — ensure BiocManager is available first
  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }

  # GitHub — ensure remotes is available first
  if (length(missing_gh) > 0) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
    for (pkg in missing_gh) {
      message("Installing from GitHub: ", packages[[pkg]]$repo)
      remotes::install_github(packages[[pkg]]$repo, upgrade = "never", quiet = TRUE)
    }
  }

  # ── Verify all packages now installed ────────────────────────────────────────

  installed_after <- rownames(installed.packages())
  still_missing   <- pkg_names[!pkg_names %in% installed_after]

  if (length(still_missing) > 0) {
    stop("The following packages could not be installed:\n",
         paste(still_missing, collapse = "\n"),
         "\nPlease install them manually and re-run the pipeline.")
  }

  # ── Load all packages ─────────────────────────────────────────────────────────

  message("Loading packages...")
  for (pkg in pkg_names) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  message("All packages loaded successfully.\n")
  invisible(TRUE)
}
