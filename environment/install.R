#!/usr/bin/env Rscript
# environment/install.R
#
# Installs the dependencies for this pipeline. Run once:
#   Rscript environment/install.R

required_cran <- c(
  "Seurat",       # core scRNA-seq analysis (v5)
  "tidyverse",    # dplyr/ggplot2/readr etc.
  "patchwork",    # combining ggplot panels
  "Matrix",       # sparse matrix handling
  "here",         # robust relative paths
  "scCustomize"   # QC + Seurat plotting helpers
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_cran, installed)

if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required CRAN packages already installed.")
}

# Record exact versions used for reproducibility
dir.create("environment", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "environment/sessionInfo.txt")
message("Session info written to environment/sessionInfo.txt")
