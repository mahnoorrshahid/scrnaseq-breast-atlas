#!/usr/bin/env Rscript
# scripts/05_compositional_analysis.R
#
# Phase 5: compare cell-type composition across the three molecular
# subtypes (ER+, HER2+, TNBC) — the piece that turns this from "I re-ran a
# tutorial" into a small original analysis with a biological conclusion.

library(Seurat)
library(tidyverse)
library(here)

obj <- readRDS(here("data", "seurat_annotated.rds"))
stopifnot(all(c("sample", "subtype", "my_celltype") %in% colnames(obj@meta.data)))

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

# ---- 1. Drop clusters excluded during annotation (low quality / doublets) -------
meta <- obj@meta.data %>%
  filter(!str_starts(my_celltype, "Excluded"))

# ---- 2. Check sample-to-subtype mapping ------------------------------------------
# Each sample (tumour) should belong to exactly one subtype.
sample_subtype <- meta %>% distinct(sample, subtype)
stopifnot(!any(duplicated(sample_subtype$sample)))
print(count(sample_subtype, subtype, name = "n_samples"))

# ---- 3. Cell-type fractions per sample --------------------------------------------
# complete() adds a 0 row when a cell type is absent from a sample;
# without it those samples would silently drop out of the comparison.
composition <- meta %>%
  count(sample, my_celltype) %>%
  complete(sample, my_celltype, fill = list(n = 0)) %>%
  group_by(sample) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  left_join(sample_subtype, by = "sample")

write_csv(composition, here("results", "05_celltype_composition_by_sample.csv"))

# ---- 4. Figures --------------------------------------------------------------------
# (a) Stacked bar: the full make-up of every tumour, grouped by subtype
p_stack <- ggplot(composition, aes(x = sample, y = fraction, fill = my_celltype)) +
  geom_col() +
  facet_grid(~subtype, scales = "free_x", space = "free_x") +
  labs(title = "Cell-type composition of each tumour",
       x = NULL, y = "Fraction of cells", fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here("figures", "05_composition_stacked_by_sample.png"), p_stack, width = 12, height = 6, dpi = 300)

# (b) Boxplots: each cell type compared across subtypes (one point per tumour)
p_box <- ggplot(composition, aes(x = subtype, y = fraction)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  facet_wrap(~my_celltype, scales = "free_y") +
  labs(title = "Cell-type composition by molecular subtype",
       x = NULL, y = "Fraction of cells per tumour") +
  theme_minimal()
ggsave(here("figures", "05_composition_by_subtype.png"), p_box, width = 12, height = 9, dpi = 300)

# ---- 5. Statistics -------------------------------------------------------------------
# Kruskal-Wallis per cell type (non-parametric: few, unevenly sized groups).
# p-values adjusted with Benjamini-Hochberg because ~10 cell types are tested.
# Caveat: fractions within a sample sum to 1, so they are not fully independent;
# a compositional method (e.g. speckle::propeller) is a possible follow-up.
test_results <- composition %>%
  group_by(my_celltype) %>%
  summarise(
    kruskal_p = tryCatch(kruskal.test(fraction ~ subtype)$p.value,
                         error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(kruskal_p, method = "BH")) %>%
  arrange(kruskal_p)

print(test_results)
write_csv(test_results, here("results", "05_composition_kruskal_tests.csv"))

# ---- 6. Exploratory: composition of the immune compartment only -------------------
# Epithelial capture varies a lot between tumours (several have ~0 epithelial
# cells, likely a dissociation artefact), which distorts whole-sample fractions.
# Re-expressing immune cell types as a fraction of immune cells removes that noise.
# Exploratory: report alongside, not instead of, the main analysis above.
immune_types <- c("T cells", "B cells", "Plasma cells", "Myeloid")

immune_comp <- composition %>%
  filter(my_celltype %in% immune_types) %>%
  group_by(sample) %>%
  mutate(frac_immune = n / sum(n)) %>%
  ungroup()

write_csv(immune_comp, here("results", "05_immune_composition_by_sample.csv"))

immune_tests <- immune_comp %>%
  group_by(my_celltype) %>%
  summarise(kruskal_p = kruskal.test(frac_immune ~ subtype)$p.value, .groups = "drop") %>%
  mutate(p_adj_BH = p.adjust(kruskal_p, method = "BH")) %>%
  arrange(kruskal_p)

print(immune_tests)
write_csv(immune_tests, here("results", "05_immune_composition_kruskal_tests.csv"))

p_immune <- ggplot(immune_comp, aes(x = subtype, y = frac_immune)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  facet_wrap(~my_celltype, scales = "free_y") +
  labs(title = "Immune compartment composition by molecular subtype",
       x = NULL, y = "Fraction of immune cells per tumour") +
  theme_minimal()
ggsave(here("figures", "05_immune_composition_by_subtype.png"), p_immune, width = 9, height = 7, dpi = 300)

# Results (n = 26 tumours: 11 ER+ / 5 HER2+ / 10 TNBC):
# - Whole-sample: no cell type significant after BH correction; myeloid nominal
#   p = 0.024 (adj 0.22), trending ER+ < HER2+ < TNBC.
# - Immune compartment: T cells nominal p = 0.019 (adj 0.077). HER2+ consistently
#   T-cell-dominated (~85% of immune cells); TNBC splits into T-cell-rich and
#   T-cell-poor tumours, the latter enriched for myeloid/plasma cells.

message("Done. Interpret cell types with p_adj_BH < 0.05 in the README's Results section.")