#!/usr/bin/env Rscript
# scripts/01_qc_filtering.R
#
# Phase 1: load the GSE176078 count matrix (26 tumours), compute per-cell QC
# metrics, filter low-quality cells, remove doublets, and save a filtered
# Seurat object for Phase 2 (integration).
#
# Required packages: Seurat, tidyverse, Matrix, here,
# scDblFinder and SingleCellExperiment (Bioconductor).

library(Seurat)
library(tidyverse)
library(Matrix)
library(here)

set.seed(42)

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

# ---- 1. Load data ----------------------------------------------------------
counts   <- readMM(here("data", "raw", "count_matrix_sparse.mtx"))
genes    <- read_tsv(here("data", "raw", "count_matrix_genes.tsv"), col_names = FALSE)$X1
barcodes <- read_tsv(here("data", "raw", "count_matrix_barcodes.tsv"), col_names = FALSE)$X1
metadata <- read_csv(here("data", "raw", "metadata.csv"))

rownames(counts) <- genes
colnames(counts) <- barcodes

message(sprintf("Loaded %d genes x %d cells (author metadata: %d cells)",
                nrow(counts), ncol(counts), nrow(metadata)))

# ---- 2. Build Seurat object -------------------------------------------------
seu <- CreateSeuratObject(
  counts = counts,
  project = "GSE176078_BrCa",
  min.cells = 3,     # gene must appear in >= 3 cells to be kept
  min.features = 200 # cell must express >= 200 genes to be kept
)

# Sample ID is the barcode prefix, e.g. "CID3838_AAACCTGAGCTAACAA-1" -> "CID3838"
seu$sample <- str_extract(colnames(seu), "^[^_]+")

# Attach author metadata (subtype, published cell-type labels) for Phases 4-5;
# not used for filtering
seu <- AddMetaData(seu, metadata = metadata %>% column_to_rownames("...1"))

# ---- 3. QC metrics ----------------------------------------------------------
seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

qc_plot <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
  group.by = "sample",
  ncol = 1,
  pt.size = 0
)
ggsave(here("figures", "01_qc_violin_by_sample.png"), qc_plot, width = 14, height = 10, dpi = 300)

# Per-sample summary, to check how much QC metrics vary between tumours
qc_summary <- seu@meta.data %>%
  group_by(sample) %>%
  summarise(
    median_nFeature = median(nFeature_RNA),
    mad_nFeature = mad(nFeature_RNA),
    median_pct_mt = median(percent_mt),
    n_cells_before = n()
  )
write_csv(qc_summary, here("results", "01_qc_summary_by_sample.csv"))

# ---- 4. Filtering ------------------------------------------------------------
# One global set of thresholds applied to all samples:
# 200-6,000 genes per cell (low = empty/broken droplet, high = likely doublet)
# and < 20% mitochondrial reads (damaged/dying cells).
# Note: no cells exceed 20% mitochondrial reads in the violin plot, suggesting
# the GEO matrix was already quality-filtered by the original authors.
seu_filtered <- subset(
  seu,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mt < 20
)

message(sprintf(
  "Cells before filtering: %d | after filtering: %d (%.1f%% retained)",
  ncol(seu), ncol(seu_filtered), 100 * ncol(seu_filtered) / ncol(seu)
))

# ---- 5. Doublet removal -------------------------------------------------------
# scDblFinder, run per sample: doublets only form between cells captured in
# the same run, so doublet rates and simulated doublets must be per sample.
library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seu_filtered)
sce <- scDblFinder(sce, samples = "sample")
seu_filtered$doublet_class <- sce$scDblFinder.class

message(sprintf("Doublets removed: %d", sum(seu_filtered$doublet_class == "doublet")))
seu_filtered <- subset(seu_filtered, subset = doublet_class == "singlet")

# ---- 6. Save -----------------------------------------------------------------
message(sprintf("Final cell count: %d across %d samples",
                ncol(seu_filtered), length(unique(seu_filtered$sample))))

saveRDS(seu_filtered, here("data", "seurat_qc_filtered.rds"))
message("Saved filtered object to data/seurat_qc_filtered.rds")
