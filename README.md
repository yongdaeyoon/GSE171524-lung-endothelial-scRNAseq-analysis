# GSE171524 Lung Endothelial scRNA-seq Analysis

This repository contains a reproducible single-cell RNA-seq analysis workflow using the public GEO dataset **GSE171524**.  
The analysis focuses on lung endothelial cell populations in control and COVID-19 lung samples, including lineage annotation, endothelial subclustering, endothelial subtype annotation, VEGFR/NRP/PTP expression, and endothelial activation signatures.

## Project goals

1. Process raw count matrices from GSE171524 into Seurat objects.
2. Perform quality control, normalization, PCA, UMAP, clustering, and marker discovery.
3. Annotate major lung lineages.
4. Subset endothelial cells and perform endothelial-only reclustering.
5. Annotate endothelial subtypes using published lung endothelial marker panels.
6. Visualize VEGFR, neuropilin, phosphatase, and endothelial activation gene expression across endothelial subtypes and disease groups.

## Dataset

- Public dataset: **GSE171524**
- Data source: NCBI Gene Expression Omnibus
- This repository does **not** redistribute raw data. Users should download raw data directly from GEO.

## Main biological focus

This workflow is designed to evaluate endothelial heterogeneity and disease-associated endothelial states in lung single-cell RNA-seq data, with emphasis on:

- Lung endothelial subtype annotation
- Arterial, venous, capillary, aerocyte-like, and lymphatic endothelial populations
- VEGFR signaling-related genes: `FLT1`, `KDR`, `FLT4`, `NRP1`, `NRP2`
- Endothelial phosphatases: `PTPRJ`, `PTPRB`, `PTPN1`
- Endothelial activation, inflammation, barrier regulation, and coagulation-related genes

## Repository structure

```text
GSE171524-lung-endothelial-scRNAseq-analysis/
├── README.md
├── .gitignore
├── R-packages.R
├── scripts/
│   ├── 01_load_qc_merge.R
│   ├── 02_whole_lung_annotation.R
│   ├── 03_endothelial_reclustering_annotation.R
│   └── 04_vegfr_activation_visualization.R
├── figures/
├── results/
└── docs/
```

## Required R packages

```r
install.packages(c("Seurat", "data.table", "Matrix", "dplyr", "ggplot2", "patchwork"))
```

## How to run

1. Download raw count files from GEO.
2. Update `raw_dir` in `scripts/01_load_qc_merge.R`.
3. Run scripts in order:

```r
source("scripts/01_load_qc_merge.R")
source("scripts/02_whole_lung_annotation.R")
source("scripts/03_endothelial_reclustering_annotation.R")
source("scripts/04_vegfr_activation_visualization.R")
```



## Endothelial-only RPCA integration workflow

This repository also includes an endothelial-only RPCA integration workflow in `scripts/05_endothelial_RPCA_integration.R`.

After whole-lung lineage annotation, endothelial cells are subsetted, split by sample, normalized independently, integrated using reciprocal PCA, reclustered, and annotated into endothelial subtypes including lymphatic, aerocyte-like, general capillary, venous, and arterial endothelial populations.

This workflow is useful for evaluating endothelial heterogeneity after reducing sample-level batch effects and for visualizing VEGFR/NRP/PTP expression across endothelial subtypes and Control/COVID groups.

## Outputs

The workflow generates:

- Whole-lung Seurat object
- Whole-lung UMAP by cluster and lineage
- Cluster marker tables
- Endothelial-only Seurat object
- Endothelial subtype UMAP
- Endothelial subtype marker DotPlots
- VEGFR/NRP/PTP expression plots
- Endothelial activation signature plots

## Skills demonstrated

- Single-cell RNA-seq analysis
- Seurat workflow
- Quality control and filtering
- Clustering and UMAP visualization
- Marker-based cell annotation
- Endothelial subtype annotation
- Differential marker analysis
- Gene module visualization
- Disease group comparison
- Reproducible bioinformatics workflow organization

## Notes

This repository is intended as a portfolio-style public-data reanalysis workflow.  
Do not upload controlled-access, unpublished, identifiable, or collaborator-owned datasets to this repository.
