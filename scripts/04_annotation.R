#!/usr/bin/env Rscript
# scripts/04_annotation.R
#
# Phase 4: annotate clusters by canonical marker genes, then cross-check
# against the original study's author-assigned labels (celltype_major,
# attached as metadata in 01_qc_filtering.R).
#
# Result: marker-based labels agreed with published labels for 92-100% of
# cells across all major lineages.

library(Seurat)
library(tidyverse)
library(here)

obj <- readRDS(here("data", "seurat_clustered.rds"))
Idents(obj) <- "res_0.8"   # clusters 0-28 (resolution 0.8) used for annotation

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

# ---- 1. Marker dot plots --------------------------------------------------------
# One canonical marker per major lineage (each gene listed once only).
major_markers <- c("EPCAM", "PTPRC", "CD3E", "MS4A1", "CD68", "PECAM1", "COL1A1", "RGS5")

p_major <- DotPlot(obj, features = major_markers) + RotatedAxis()
ggsave(here("figures", "04_dotplot_major_markers.png"), p_major, width = 9, height = 8, dpi = 300)

# Follow-up panel for clusters that were unclear from the major panel.
followup_markers <- c("MZB1", "JCHAIN", "PTPRC", "PECAM1", "RGS5", "ACTA2",
                      "EPCAM", "KRT19", "MKI67", "NKG7")
unclear_clusters <- c(11, 15, 21, 23, 24, 25, 27)

p_followup <- DotPlot(obj, features = followup_markers, idents = unclear_clusters) + RotatedAxis()
ggsave(here("figures", "04_dotplot_unclear_clusters.png"), p_followup, width = 9, height = 6, dpi = 300)

# ---- 2. Marker-based labels (before validation) ----------------------------------
marker_labels <- c(
  "0"  = "Epithelial",   "1"  = "T cells",      "2"  = "T cells",
  "3"  = "T cells",      "4"  = "CAFs",         "5"  = "Endothelial",
  "6"  = "Myeloid",      "7"  = "PVL",          "8"  = "Epithelial",
  "9"  = "Epithelial",   "10" = "Epithelial",   "11" = "Plasma cells",   # MZB1/JCHAIN+
  "12" = "B cells",      "13" = "Myeloid",      "14" = "T cells",
  "15" = "Unclear (ACTA2+)",                                             # ACTA2+ only
  "16" = "Endothelial",  "17" = "Epithelial",   "18" = "T cells",
  "19" = "Epithelial",   "20" = "Myeloid",
  "21" = "Unclear (vascular)",                                           # PECAM1+ RGS5+
  "22" = "Myeloid",
  "23" = "Low quality?",                                                 # weak for all markers
  "24" = "Plasma cells",                                                 # MZB1/JCHAIN/MKI67+ (plasmablasts)
  "25" = "Low quality?",                                                 # weak for all markers
  "26" = "Epithelial",
  "27" = "Doublets?",                                                    # EPCAM/KRT19 + RGS5/ACTA2
  "28" = "Epithelial"
)
obj$marker_celltype <- unname(marker_labels[as.character(Idents(obj))])

# ---- 3. Validate against published labels ----------------------------------------
agreement <- round(prop.table(table(mine = obj$marker_celltype,
                                    published = obj$celltype_major), 1) * 100)
print(agreement)
write.csv(as.data.frame.matrix(agreement), here("results", "04_annotation_vs_published_pct.csv"))

# ---- 4. Final labels (unclear clusters resolved using the validation table) ------
# Unclear (ACTA2+)   -> 94% Normal Epithelial in published labels: myoepithelial cells
# Unclear (vascular) -> 98% Endothelial
# Low quality?       -> mixed published identity (44% cancer epi / 34% endo / 17% T): excluded
# Doublets?          -> co-expresses epithelial and PVL markers: excluded as likely doublets
obj$my_celltype <- recode(obj$marker_celltype,
                          "Unclear (ACTA2+)"   = "Normal epithelial (myoepithelial)",
                          "Unclear (vascular)" = "Endothelial",
                          "Low quality?"       = "Excluded (low quality)",
                          "Doublets?"          = "Excluded (doublets)"
)
print(table(obj$my_celltype))

p_annotated <- DimPlot(obj, reduction = "umap.integrated", group.by = "my_celltype", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Annotated cell types")
ggsave(here("figures", "04_umap_annotated.png"), p_annotated, width = 9, height = 7, dpi = 300)

# ---- 5. Save -----------------------------------------------------------------------
saveRDS(obj, here("data", "seurat_annotated.rds"))
message("Saved annotated object to data/seurat_annotated.rds")
