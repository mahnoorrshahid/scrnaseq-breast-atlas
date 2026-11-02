#!/usr/bin/env Rscript
# scripts/03_clustering.R
#
# Phase 3: graph-based clustering on the Harmony-integrated embedding,
# choice of resolution with clustree, and marker genes for each cluster.
#
# Required packages: Seurat, clustree, tidyverse, here.

library(Seurat)
library(clustree)
library(tidyverse)
library(here)

set.seed(42)

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

seu <- readRDS(here("data", "seurat_integrated.rds"))

# ---- 1. Clustering -----------------------------------------------------------
seu <- FindNeighbors(seu, reduction = "integrated.harmony", dims = 1:30)

# Several resolutions compared rather than one assumed: too low merges distinct
# populations (e.g. CAFs vs PVL cells), too high fragments one cell type into
# near-duplicate clusters.
resolutions <- c(0.4, 0.6, 0.8, 1.0)
seu <- FindClusters(seu, resolution = resolutions, cluster.name = paste0("res_", resolutions))

# Clustering tree: how clusters at each resolution split at the next one
ct <- clustree(seu@meta.data, prefix = "res_")
ggsave(here("figures", "03_clustree.png"), ct, width = 10, height = 12, dpi = 300)

# Resolution 0.8 chosen: major clusters stay stable from 0.4 to 0.8, while at
# 1.0 several large clusters fragment. Gives 29 clusters (0-28).
Idents(seu) <- "res_0.8"

seu <- RunUMAP(seu, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.clustered")

p_clusters <- DimPlot(seu, reduction = "umap.clustered", label = TRUE) +
  ggtitle("Clusters (resolution 0.8)")
ggsave(here("figures", "03_umap_clusters.png"), p_clusters, width = 8, height = 7, dpi = 300)

# ---- 2. Marker genes ----------------------------------------------------------
# Genes up-regulated in each cluster vs all other cells: positive markers only,
# expressed in >= 25% of the cluster's cells, log2 fold-change >= 0.25
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write_csv(markers, here("results", "03_all_cluster_markers.csv"))
write_csv(top_markers, here("results", "03_top10_markers_per_cluster.csv"))

marker_heatmap <- DoHeatmap(seu, features = unique(top_markers$gene)) + NoLegend()
ggsave(here("figures", "03_marker_heatmap.png"), marker_heatmap, width = 24, height = 20, dpi = 300)

# ---- 3. Save -------------------------------------------------------------------
saveRDS(seu, here("data", "seurat_clustered.rds"))
message(sprintf("Identified %d clusters at resolution 0.8.", length(unique(seu$res_0.8))))
message("Saved clustered object to data/seurat_clustered.rds")
