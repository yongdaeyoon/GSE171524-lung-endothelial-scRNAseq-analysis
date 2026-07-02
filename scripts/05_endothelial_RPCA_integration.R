# 05_endothelial_RPCA_integration.R
# Endothelial-only RPCA integration workflow for GSE171524.
#
# Purpose:
# After whole-lung lineage annotation, endothelial cells are subsetted,
# split by sample, integrated using reciprocal PCA (RPCA), reclustered,
# annotated by endothelial subtype markers, and visualized for VEGFR/NRP/PTP
# expression across endothelial subtypes and Control/COVID groups.
#
# Required input:
# results/GSE171524_whole_lung_annotated.rds
#
# Expected output:
# results/GSE171524_endothelial_RPCA_integrated_annotated.rds
# figures/endothelial_RPCA_*.pdf

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------
# Paths
# -----------------------------
output_dir <- "results"
figure_dir <- "figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load whole-lung annotated object
# -----------------------------
covid <- readRDS(file.path(output_dir, "GSE171524_whole_lung_annotated.rds"))

DefaultAssay(covid) <- "RNA"

# -----------------------------
# Subset endothelial cells
# -----------------------------
endo <- subset(
  covid,
  subset = lineage_manual == "Endothelial"
)

message("Initial endothelial subset:")
print(endo)
print(table(endo$group))
print(table(endo$sample_id))

# -----------------------------
# Optional stringent endothelial marker-positive filtering
# -----------------------------
endo <- subset(
  endo,
  subset =
    PECAM1 > 0 |
    CDH5 > 0 |
    CLDN5 > 0 |
    VWF > 0 |
    EMCN > 0 |
    GJA5 > 0 |
    BMX > 0 |
    HEY1 > 0 |
    SEMA3G > 0
)

message("Endothelial marker-positive subset:")
print(endo)
print(table(endo$group))
print(table(endo$sample_id))

saveRDS(
  endo,
  file = file.path(output_dir, "GSE171524_endothelial_before_RPCA_integration.rds")
)

# -----------------------------
# Split endothelial object by sample
# -----------------------------
DefaultAssay(endo) <- "RNA"

endo_list <- SplitObject(
  endo,
  split.by = "sample_id"
)

message("Number of endothelial sample objects:")
print(length(endo_list))

message("Cells per endothelial sample object:")
print(sapply(endo_list, ncol))

# -----------------------------
# Normalize and find variable features per sample
# -----------------------------
endo_list <- lapply(endo_list, function(x) {
  DefaultAssay(x) <- "RNA"

  x <- NormalizeData(x)

  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )

  return(x)
})

# -----------------------------
# Select integration features
# -----------------------------
features_endo <- SelectIntegrationFeatures(
  object.list = endo_list,
  nfeatures = 2000
)

# -----------------------------
# Prepare each object for RPCA integration
# -----------------------------
endo_list <- lapply(endo_list, function(x) {
  x <- ScaleData(
    x,
    features = features_endo,
    verbose = FALSE
  )

  x <- RunPCA(
    x,
    features = features_endo,
    npcs = 30,
    verbose = FALSE
  )

  return(x)
})

# -----------------------------
# Find integration anchors using RPCA
# -----------------------------
anchors_endo <- FindIntegrationAnchors(
  object.list = endo_list,
  anchor.features = features_endo,
  reduction = "rpca",
  dims = 1:20
)

saveRDS(
  anchors_endo,
  file = file.path(output_dir, "GSE171524_endothelial_RPCA_anchors.rds")
)

# -----------------------------
# Integrate endothelial cells
# -----------------------------
endo_integrated <- IntegrateData(
  anchorset = anchors_endo,
  dims = 1:20,
  k.weight = 30
)

# -----------------------------
# Integrated PCA / UMAP / clustering
# -----------------------------
DefaultAssay(endo_integrated) <- "integrated"

endo_integrated <- ScaleData(
  endo_integrated,
  verbose = FALSE
)

endo_integrated <- RunPCA(
  endo_integrated,
  npcs = 50,
  verbose = FALSE
)

pdf(file.path(figure_dir, "endothelial_RPCA_elbowplot.pdf"), width = 6, height = 4)
print(ElbowPlot(endo_integrated, ndims = 50))
dev.off()

endo_integrated <- FindNeighbors(
  endo_integrated,
  dims = 1:15
)

endo_integrated <- FindClusters(
  endo_integrated,
  resolution = 0.3
)

endo_integrated <- RunUMAP(
  endo_integrated,
  dims = 1:15
)

# -----------------------------
# Integration QC plots
# -----------------------------
p_clusters <- DimPlot(
  endo_integrated,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("EC-only integrated clusters")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_integrated_clusters_umap.pdf"),
  p_clusters,
  width = 7,
  height = 5
)

p_sample <- DimPlot(
  endo_integrated,
  reduction = "umap",
  group.by = "sample_id",
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("EC-only integration by sample")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_integrated_by_sample_umap.pdf"),
  p_sample,
  width = 8,
  height = 5
)

p_group <- DimPlot(
  endo_integrated,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("Control vs COVID in EC-only integrated UMAP")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_integrated_by_group_umap.pdf"),
  p_group,
  width = 7,
  height = 5
)

# -----------------------------
# Endothelial subtype marker panels
# Marker expression should be checked on RNA assay, not integrated assay.
# -----------------------------
DefaultAssay(endo_integrated) <- "RNA"
Idents(endo_integrated) <- "seurat_clusters"

ec_check_markers <- list(
  Lymphatic = c("PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "RELN"),
  Aerocyte = c("CA4", "EDNRB", "SOSTDC1", "HPGD", "TBX2", "CYP3A5"),
  gCap = c("RGCC", "FCN3", "GPIHBP1", "CD36", "BTNL9", "ADGRF5"),
  Vein = c("ACKR1", "VWF", "PLVAP", "COL15A1", "CLU", "C7", "PLAT"),
  Arterial = c("GJA5", "BMX", "HEY1", "SEMA3G", "DKK2", "CXCL12")
)

genes_use <- unique(unlist(ec_check_markers))
genes_use <- genes_use[genes_use %in% rownames(endo_integrated)]

p_marker_clusters <- DotPlot(
  endo_integrated,
  features = genes_use,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("EC subtype marker check after RPCA integration")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_subtype_marker_by_cluster_dotplot.pdf"),
  p_marker_clusters,
  width = 12,
  height = 5
)

# -----------------------------
# Manual endothelial subtype annotation after integration
# Update these cluster IDs after inspecting marker plots.
# -----------------------------
ec_lymph_clusters    <- c("1", "8")
ec_aerocyte_clusters <- c("0")
ec_gcap_clusters     <- c("3", "5", "6", "7", "9")
ec_vein_clusters     <- c("4")
ec_arterial_clusters <- c("2")

endo_integrated$ec_subtype_manual <- "Unknown"

endo_integrated$ec_subtype_manual[
  as.character(Idents(endo_integrated)) %in% ec_lymph_clusters
] <- "Lymphatic"

endo_integrated$ec_subtype_manual[
  as.character(Idents(endo_integrated)) %in% ec_aerocyte_clusters
] <- "Aerocyte"

endo_integrated$ec_subtype_manual[
  as.character(Idents(endo_integrated)) %in% ec_gcap_clusters
] <- "gCap"

endo_integrated$ec_subtype_manual[
  as.character(Idents(endo_integrated)) %in% ec_vein_clusters
] <- "Vein"

endo_integrated$ec_subtype_manual[
  as.character(Idents(endo_integrated)) %in% ec_arterial_clusters
] <- "Arterial"

message("Endothelial subtype counts after RPCA integration:")
print(table(endo_integrated$ec_subtype_manual))

p_subtype <- DimPlot(
  endo_integrated,
  reduction = "umap",
  group.by = "ec_subtype_manual",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("EC subtype annotation after RPCA integration")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_subtype_annotation_umap.pdf"),
  p_subtype,
  width = 7,
  height = 5
)

p_marker_subtype <- DotPlot(
  endo_integrated,
  features = genes_use,
  group.by = "ec_subtype_manual",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("EC subtype marker check by manual annotation")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_subtype_marker_by_annotation_dotplot.pdf"),
  p_marker_subtype,
  width = 12,
  height = 5
)

# -----------------------------
# Subtype marker FeaturePlots
# -----------------------------
feature_groups <- list(
  Lymphatic = c("PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "RELN"),
  Aerocyte = c("CA4", "EDNRB", "SOSTDC1", "HPGD", "TBX2", "CYP3A5"),
  gCap = c("RGCC", "FCN3", "GPIHBP1", "CD36", "BTNL9", "ADGRF5"),
  Vein = c("ACKR1", "VWF", "PLVAP", "COL15A1", "CLU", "C7", "PLAT"),
  Arterial = c("GJA5", "BMX", "HEY1", "SEMA3G", "DKK2", "CXCL12")
)

for (nm in names(feature_groups)) {
  genes <- intersect(feature_groups[[nm]], rownames(endo_integrated))
  if (length(genes) > 0) {
    p <- FeaturePlot(
      endo_integrated,
      features = genes,
      reduction = "umap",
      ncol = 3,
      pt.size = 0.4,
      order = TRUE
    ) +
      ggtitle(paste0(nm, " marker genes"))

    ggsave(
      file.path(figure_dir, paste0("endothelial_RPCA_featureplot_", nm, "_markers.pdf")),
      p,
      width = 10,
      height = 6
    )
  }
}

# -----------------------------
# VEGFR / NRP / PTP expression after RPCA integration
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

vegfr_ptp_genes_present <- vegfr_ptp_genes[
  vegfr_ptp_genes %in% rownames(endo_integrated)
]

vegfr_ptp_genes_missing <- setdiff(
  vegfr_ptp_genes,
  rownames(endo_integrated)
)

message("VEGFR/NRP/PTP genes present:")
print(vegfr_ptp_genes_present)

message("VEGFR/NRP/PTP genes missing:")
print(vegfr_ptp_genes_missing)

p_vegfr_subtype <- DotPlot(
  endo_integrated,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_manual",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("VEGFR / NRP / PTP expression by EC subtype after RPCA integration")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_vegfr_nrp_ptp_by_subtype_dotplot.pdf"),
  p_vegfr_subtype,
  width = 8,
  height = 4
)

endo_integrated$ec_subtype_group <- paste(
  endo_integrated$ec_subtype_manual,
  endo_integrated$group,
  sep = "_"
)

p_vegfr_subtype_group <- DotPlot(
  endo_integrated,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_group",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("VEGFR / NRP / PTP expression by EC subtype and group")

ggsave(
  file.path(figure_dir, "endothelial_RPCA_vegfr_nrp_ptp_by_subtype_group_dotplot.pdf"),
  p_vegfr_subtype_group,
  width = 10,
  height = 5
)

p_vegfr_feature <- FeaturePlot(
  endo_integrated,
  features = vegfr_ptp_genes_present,
  reduction = "umap",
  ncol = 4,
  pt.size = 0.4,
  order = TRUE
)

ggsave(
  file.path(figure_dir, "endothelial_RPCA_vegfr_nrp_ptp_featureplots.pdf"),
  p_vegfr_feature,
  width = 12,
  height = 6
)

p_vegfr_vln <- VlnPlot(
  endo_integrated,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_manual",
  pt.size = 0.1,
  ncol = 3
) &
  theme_classic()

ggsave(
  file.path(figure_dir, "endothelial_RPCA_vegfr_nrp_ptp_vlnplot_by_subtype.pdf"),
  p_vegfr_vln,
  width = 10,
  height = 7
)

p_vegfr_vln_split <- VlnPlot(
  endo_integrated,
  features = vegfr_ptp_genes_present,
  group.by = "ec_subtype_manual",
  split.by = "group",
  pt.size = 0.1,
  ncol = 3
) &
  theme_classic()

ggsave(
  file.path(figure_dir, "endothelial_RPCA_vegfr_nrp_ptp_vlnplot_by_subtype_group.pdf"),
  p_vegfr_vln_split,
  width = 10,
  height = 7
)

# -----------------------------
# Save final integrated object
# -----------------------------
saveRDS(
  endo_integrated,
  file = file.path(output_dir, "GSE171524_endothelial_RPCA_integrated_annotated.rds")
)

message("Done: endothelial RPCA integration workflow completed.")
