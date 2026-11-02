#!/usr/bin/env Rscript
# scripts/02_integration.R
#
# Phase 2: normalise each sample, find variable genes, run PCA, and integrate
# across the 26 patients with Harmony to correct for patient-driven batch
# effects before clustering.
#
# Required packages: Seurat (v5), harmony, tidyverse, patchwork, here.

library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(here)

set.seed(42)

dir.create(here("figures"), showWarnings = FALSE)

seu <- readRDS(here("data", "seurat_qc_filtered.rds"))

# ---- 1. Split by sample, normalise, PCA -------------------------------------
# Seurat v5: splitting the RNA assay into one layer per sample lets
# IntegrateLayers() treat each tumour as a separate batch.
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)

seu <- NormalizeData(seu)                       # log-normalisation (default)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)

# ---- 2. Integrate with Harmony ------------------------------------------------
# Seurat's CCA integration was tried first but was too slow to complete on
# ~98,000 cells, so Harmony was used: it is designed for large datasets and
# corrects the PCA embedding without altering the underlying counts.
seu <- IntegrateLayers(
  object = seu,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  verbose = TRUE
)

# Re-join the per-sample layers for downstream marker finding
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

# ---- 3. Check: UMAP before vs after integration -------------------------------
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
seu <- RunUMAP(seu, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.integrated")

p1 <- DimPlot(seu, reduction = "umap.unintegrated", group.by = "sample") +
  ggtitle("Before integration (by sample)") + NoLegend()
p2 <- DimPlot(seu, reduction = "umap.integrated", group.by = "sample") +
  ggtitle("After integration (by sample)") + NoLegend()

ggsave(here("figures", "02_integration_before_after.png"),
       p1 + p2, width = 14, height = 6, dpi = 300)

# ---- 4. Save -------------------------------------------------------------------
saveRDS(seu, here("data", "seurat_integrated.rds"))
message("Saved integrated object to data/seurat_integrated.rds")
