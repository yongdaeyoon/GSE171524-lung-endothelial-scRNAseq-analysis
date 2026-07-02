# 02_whole_lung_annotation.R
# Whole-lung normalization, clustering, marker detection, and major lineage annotation.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

output_dir <- "results"
figure_dir <- "figures"

covid <- readRDS(file.path(output_dir, "GSE171524_merged_QC.rds"))

DefaultAssay(covid) <- "RNA"

# -----------------------------
# Normalize, PCA, UMAP, clustering
# -----------------------------
covid <- NormalizeData(covid)
covid <- FindVariableFeatures(covid, selection.method = "vst", nfeatures = 3000)
covid <- ScaleData(covid, features = VariableFeatures(covid))
covid <- RunPCA(covid, features = VariableFeatures(covid), npcs = 50)

pdf(file.path(figure_dir, "whole_lung_elbowplot.pdf"), width = 6, height = 4)
print(ElbowPlot(covid, ndims = 50))
dev.off()

covid <- FindNeighbors(covid, dims = 1:30)
covid <- FindClusters(covid, resolution = 0.5)
covid <- RunUMAP(covid, dims = 1:30)

p_cluster <- DimPlot(
  covid,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  theme_classic() +
  ggtitle("GSE171524 whole-lung Seurat clusters")

ggsave(file.path(figure_dir, "whole_lung_clusters_umap.pdf"), p_cluster, width = 7, height = 5)

# -----------------------------
# Cluster markers
# -----------------------------
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  covid <- JoinLayers(covid)
}

Idents(covid) <- "seurat_clusters"

markers_all <- FindAllMarkers(
  covid,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10
)

write.csv(
  markers_all,
  file = file.path(output_dir, "GSE171524_whole_lung_cluster_markers.csv"),
  row.names = FALSE
)

top10_all <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top10_all,
  file = file.path(output_dir, "GSE171524_whole_lung_cluster_top10_markers.csv"),
  row.names = FALSE
)

# -----------------------------
# Major lineage marker check
# -----------------------------
lineage_markers <- c(
  "LYZ", "C1QA", "C1QB", "CD68",
  "PTPRC", "CD3D", "CD3E", "NKG7", "MS4A1", "CD79A",
  "EPCAM", "CDH1", "KRT18", "KRT19", "AGER", "SFTPC",
  "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB",
  "PECAM1", "CLDN5", "CDH5", "VWF", "EMCN"
)

lineage_markers <- lineage_markers[lineage_markers %in% rownames(covid)]

p_lineage_dot <- DotPlot(covid, features = lineage_markers, dot.scale = 6) +
  RotatedAxis() +
  theme_classic() +
  ggtitle("Major lineage marker check")

ggsave(file.path(figure_dir, "whole_lung_lineage_marker_dotplot.pdf"), p_lineage_dot, width = 12, height = 5)

# -----------------------------
# Manual lineage annotation
# Update cluster IDs after inspecting markers if needed.
# -----------------------------
myeloid_clusters     <- c("0", "6", "9", "20")
lymphoid_clusters    <- c("1", "8", "11", "17", "18")
epithelial_clusters  <- c("3", "4", "5", "12", "13", "15")
stromal_clusters     <- c("2", "10", "14", "19", "21", "22")
endothelial_clusters <- c("7", "16")

covid$lineage_manual <- "Unknown"

covid$lineage_manual[as.character(Idents(covid)) %in% myeloid_clusters]     <- "Myeloid"
covid$lineage_manual[as.character(Idents(covid)) %in% lymphoid_clusters]    <- "Lymphoid"
covid$lineage_manual[as.character(Idents(covid)) %in% epithelial_clusters]  <- "Epithelial"
covid$lineage_manual[as.character(Idents(covid)) %in% stromal_clusters]     <- "Stromal"
covid$lineage_manual[as.character(Idents(covid)) %in% endothelial_clusters] <- "Endothelial"

p_lineage <- DimPlot(
  covid,
  reduction = "umap",
  group.by = "lineage_manual",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  theme_classic() +
  ggtitle("GSE171524 major lineage annotation")

ggsave(file.path(figure_dir, "whole_lung_lineage_annotation_umap.pdf"), p_lineage, width = 7, height = 5)

saveRDS(covid, file = file.path(output_dir, "GSE171524_whole_lung_annotated.rds"))
