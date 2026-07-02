# 01_load_qc_merge.R
# Load raw count matrices, create Seurat objects, perform QC filtering, and merge samples.

library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------
# User settings
# -----------------------------
raw_dir <- "PATH/TO/GSE171524_RAW"   # Update this path
output_dir <- "results"
figure_dir <- "figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Find raw count files
# -----------------------------
csv_files <- list.files(
  raw_dir,
  pattern = "raw_counts\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(length(csv_files) > 0)

# -----------------------------
# Function: read one sample
# -----------------------------
read_one_sample <- function(file) {
  sample_id <- basename(file)
  sample_id <- gsub("_raw_counts.csv", "", sample_id)
  sample_id <- gsub("^GSM[0-9]+_", "", sample_id)

  message("Reading: ", sample_id)

  mat <- fread(file, data.table = FALSE)
  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]

  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat <- Matrix(mat, sparse = TRUE)

  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )

  obj$sample_id <- sample_id
  obj$group <- ifelse(grepl("ctr", sample_id, ignore.case = TRUE), "Control", "COVID")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  return(obj)
}

# -----------------------------
# Read all samples
# -----------------------------
obj_list <- lapply(csv_files, read_one_sample)

# -----------------------------
# QC filtering
# -----------------------------
obj_list <- lapply(obj_list, function(x) {
  subset(
    x,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 6000 &
      percent.mt < 20
  )
})

# -----------------------------
# Merge samples
# -----------------------------
covid <- merge(
  obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = sapply(obj_list, function(x) unique(x$sample_id)),
  project = "GSE171524"
)

message("Cells by group:")
print(table(covid$group))

message("Cells by sample:")
print(table(covid$sample_id))

saveRDS(covid, file = file.path(output_dir, "GSE171524_merged_QC.rds"))
