# 04_vegfr_activation_visualization.R
# Visualize VEGFR/NRP/PTP genes and endothelial activation signatures.

library(Seurat)
library(ggplot2)
library(patchwork)

output_dir <- "results"
figure_dir <- "figures"

endo <- readRDS(file.path(output_dir, "GSE171524_endothelial_annotated.rds"))

DefaultAssay(endo) <- "RNA"

# -----------------------------
# VEGFR / NRP / PTP gene panel
# -----------------------------
vegfr_ptp_genes <- c(
  "FLT1",   # VEGFR1
  "KDR",    # VEGFR2
  "FLT4",   # VEGFR3
  "NRP1",
  "NRP2",
  "PTPRJ",  # DEP-1
  "PTPRB",  # VE-PTP
  "PTPN1"   # PTP1B
)

vegfr_ptp_genes_present <- intersect(vegfr_ptp_genes, rownames(endo))
vegfr_ptp_genes_missing <- setdiff(vegfr_ptp_genes, rownames(endo))

writeLines("Present genes:")
print(vegfr_ptp_genes_present)

writeLines("Missing genes:")
print(vegfr_ptp_genes_missing)

p_dot_subtype <- DotPlot(
  endo,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_manual",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("VEGFR / NRP / PTP expression by EC subtype")

ggsave(file.path(figure_dir, "vegfr_nrp_ptp_by_ec_subtype_dotplot.pdf"), p_dot_subtype, width = 8, height = 4)

endo$ec_subtype_group <- paste(endo$ec_subtype_manual, endo$group, sep = "_")

p_dot_subtype_group <- DotPlot(
  endo,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_group",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("VEGFR / NRP / PTP expression by EC subtype and group")

ggsave(file.path(figure_dir, "vegfr_nrp_ptp_by_ec_subtype_group_dotplot.pdf"), p_dot_subtype_group, width = 10, height = 5)

# -----------------------------
# Helper function: gene-positive cells colored by group
# -----------------------------
plot_gene_group_umap <- function(seurat_obj, gene) {
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

  umap_df$group <- seurat_obj$group
  umap_df$expr <- FetchData(seurat_obj, vars = gene)[, 1]

  ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      color = "lightgrey",
      size = 0.2,
      alpha = 0.4
    ) +
    geom_point(
      data = subset(umap_df, expr > 0),
      aes(color = group, size = expr),
      alpha = 0.8
    ) +
    scale_color_manual(values = c("Control" = "blue", "COVID" = "red")) +
    scale_size_continuous(range = c(0.2, 1.2)) +
    theme_classic() +
    ggtitle(paste0(gene, " expression: Control vs COVID"))
}

for (gene in vegfr_ptp_genes_present) {
  p <- plot_gene_group_umap(endo, gene)
  ggsave(file.path(figure_dir, paste0("umap_", gene, "_control_vs_covid.pdf")), p, width = 6, height = 5)
}

# -----------------------------
# Endothelial activation gene panel
# -----------------------------
ec_activation <- c(
  # Rolling and adhesion
  "ICAM1", "VCAM1", "SELE", "SELP", "CXCL12", "CCL2",

  # IFN-stimulated genes
  "CXCL10", "CXCL9", "ISG15", "IFI6", "MX1", "STAT1",

  # Thrombo / coagulation
  "SERPINE1", "F3", "PLAT", "THBD",

  # Barrier dysfunction and EndoMT
  "ANGPT2", "CDH5", "CLDN5", "TJP1", "RAMP2", "ADM", "ADM2", "TEK", "CALCRL",

  # Cytokine-related
  "IL1B", "TNF", "IL1R1", "IL6", "TGFB1", "TGFB2"
)

ec_activation_present <- intersect(ec_activation, rownames(endo))
ec_activation_missing <- setdiff(ec_activation, rownames(endo))

writeLines("Activation genes present:")
print(ec_activation_present)

writeLines("Activation genes missing:")
print(ec_activation_missing)

p_activation_dot <- DotPlot(
  endo,
  features = ec_activation_present,
  group.by = "ec_subtype_manual",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("Endothelial activation genes by EC subtype")

ggsave(file.path(figure_dir, "endothelial_activation_by_subtype_dotplot.pdf"), p_activation_dot, width = 14, height = 5)

# -----------------------------
# Module score example
# -----------------------------
ec_gene_sets <- list(
  Rolling_adhesion = c("ICAM1", "VCAM1", "SELE", "SELP", "CXCL12", "CCL2"),
  IFN_response = c("CXCL10", "CXCL9", "ISG15", "IFI6", "MX1", "STAT1"),
  Thrombo = c("SERPINE1", "F3", "PLAT", "THBD"),
  Barrier_EndoMT = c("ANGPT2", "CDH5", "CLDN5", "TJP1", "RAMP2", "ADM", "ADM2", "TEK", "CALCRL"),
  Cytokine = c("IL1B", "TNF", "IL1R1", "IL6", "TGFB1", "TGFB2")
)

ec_gene_sets_present <- lapply(ec_gene_sets, function(x) intersect(x, rownames(endo)))
ec_gene_sets_present <- ec_gene_sets_present[lengths(ec_gene_sets_present) > 0]

endo <- AddModuleScore(
  object = endo,
  features = ec_gene_sets_present,
  name = names(ec_gene_sets_present)
)

# AddModuleScore appends numbers to names in order.
score_cols <- grep(
  paste(names(ec_gene_sets_present), collapse = "|"),
  colnames(endo@meta.data),
  value = TRUE
)

p_score_group <- VlnPlot(
  endo,
  features = score_cols,
  group.by = "group",
  pt.size = 0,
  ncol = 3
) &
  theme_classic()

ggsave(file.path(figure_dir, "endothelial_activation_module_scores_by_group.pdf"), p_score_group, width = 10, height = 6)

saveRDS(endo, file = file.path(output_dir, "GSE171524_endothelial_annotated_with_scores.rds"))
