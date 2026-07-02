# 03_endothelial_reclustering_annotation.R
# Subset endothelial cells, recluster, and annotate endothelial subtypes.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

output_dir <- "results"
figure_dir <- "figures"

covid <- readRDS(file.path(output_dir, "GSE171524_whole_lung_annotated.rds"))

# -----------------------------
# Subset endothelial cells
# -----------------------------
endo <- subset(covid, subset = lineage_manual == "Endothelial")

# Optional stricter endothelial marker-positive filter
endo <- subset(
  endo,
  subset =
    PECAM1 > 0 | CDH5 > 0 | CLDN5 > 0 | VWF > 0 | EMCN > 0 |
    GJA5 > 0 | BMX > 0 | HEY1 > 0 | SEMA3G > 0
)

message("Endothelial cells by group:")
print(table(endo$group))

# -----------------------------
# EC-only preprocessing
# -----------------------------
DefaultAssay(endo) <- "RNA"

endo <- NormalizeData(endo)
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 3000)
endo <- ScaleData(endo, features = VariableFeatures(endo))
endo <- RunPCA(endo, features = VariableFeatures(endo), npcs = 50)

pdf(file.path(figure_dir, "endothelial_elbowplot.pdf"), width = 6, height = 4)
print(ElbowPlot(endo, ndims = 50))
dev.off()

endo <- FindNeighbors(endo, dims = 1:15)
endo <- FindClusters(endo, resolution = 0.3)
endo <- RunUMAP(endo, dims = 1:15)

p_ec_cluster <- DimPlot(
  endo,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("Endothelial subclusters: dims 1:15, res 0.3")

ggsave(file.path(figure_dir, "endothelial_subclusters_umap.pdf"), p_ec_cluster, width = 7, height = 5)

# -----------------------------
# EC subtype marker panel
# -----------------------------
ec_markers_by_type <- list(
  EC_arterial = c(
    "DKK2", "IGFBP3", "FBLN5", "SERPINE2", "CLDN10",
    "GJA5", "CXCL12", "BMX", "LTBP4", "HEY1",
    "SOX5", "SEMA3G"
  ),

  EC_aerocyte = c(
    "SOSTDC1", "EDNRB", "HPGD", "CYP3A5", "PRKG1",
    "TBX2", "RCSD1", "EDA", "B3GALNT1", "EXPH5",
    "NCALD", "S100A4", "CA4", "AFF3", "ADGRL2"
  ),

  EC_general_capillary = c(
    "BTNL9", "RGCC", "ADGRF5", "KIAA1217", "FCN3",
    "IL7R", "CD36", "NRXN3", "SLC6A4", "GPIHBP1",
    "ARHGAP18"
  ),

  EC_venous = c(
    "ACKR1", "VWF", "PLVAP", "COL15A1", "IGFBP7",
    "CLU", "C7", "PLAT", "DKK3", "CPE", "PTGS1",
    "MMRN1", "PKHD1L1"
  ),

  EC_lymphatic = c(
    "PROX1", "LYVE1", "PDPN", "FLT4", "CCL21",
    "SEMA3D", "TFF3", "TM4SF18", "TBX1", "RELN"
  )
)

ec_markers_by_type_present <- lapply(
  ec_markers_by_type,
  function(x) intersect(x, rownames(endo))
)

p_ec_dot <- DotPlot(
  endo,
  features = unique(unlist(ec_markers_by_type_present)),
  group.by = "seurat_clusters",
  dot.scale = 5
) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("Endothelial subtype marker panel")

ggsave(file.path(figure_dir, "endothelial_subtype_marker_dotplot.pdf"), p_ec_dot, width = 14, height = 6)

# -----------------------------
# Manual EC subtype annotation
# Update cluster IDs after inspecting markers if needed.
# -----------------------------
Idents(endo) <- "seurat_clusters"

ec_capillary_clusters <- c("0", "5", "6")
ec_lymphatic_clusters <- c("1", "7")
ec_venous_clusters    <- c("2")
ec_aerocyte_clusters  <- c("3")
ec_arterial_clusters  <- c("4")

endo$ec_subtype_manual <- "Unknown / mixed EC"

endo$ec_subtype_manual[as.character(Idents(endo)) %in% ec_capillary_clusters] <- "gCap"
endo$ec_subtype_manual[as.character(Idents(endo)) %in% ec_lymphatic_clusters] <- "Lymphatic"
endo$ec_subtype_manual[as.character(Idents(endo)) %in% ec_venous_clusters]    <- "Venous"
endo$ec_subtype_manual[as.character(Idents(endo)) %in% ec_aerocyte_clusters]  <- "Aerocyte"
endo$ec_subtype_manual[as.character(Idents(endo)) %in% ec_arterial_clusters]  <- "Arterial"

message("Endothelial subtype counts:")
print(table(endo$ec_subtype_manual))

p_ec_annot <- DimPlot(
  endo,
  reduction = "umap",
  group.by = "ec_subtype_manual",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("Manual endothelial subtype annotation")

ggsave(file.path(figure_dir, "endothelial_subtype_annotation_umap.pdf"), p_ec_annot, width = 7, height = 5)

# -----------------------------
# EC cluster markers
# -----------------------------
endo_markers <- FindAllMarkers(
  endo,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10
)

write.csv(
  endo_markers,
  file = file.path(output_dir, "GSE171524_endothelial_cluster_markers.csv"),
  row.names = FALSE
)

saveRDS(endo, file = file.path(output_dir, "GSE171524_endothelial_annotated.rds"))
