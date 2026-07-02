# Endothelial-only RPCA Integration Notes

This workflow extends the main GSE171524 analysis by performing integration after endothelial cell subsetting.

Main steps:

1. Load whole-lung annotated object.
2. Subset endothelial cells.
3. Apply a stringent endothelial marker-positive filter.
4. Split endothelial cells by sample.
5. Normalize each endothelial sample object.
6. Select integration features.
7. Run PCA per sample object.
8. Find integration anchors using RPCA.
9. Integrate data with `k.weight = 30`.
10. Run PCA, UMAP, and clustering on the integrated assay.
11. Switch back to RNA assay for marker validation.
12. Annotate EC subtypes using marker panels.
13. Visualize VEGFR/NRP/PTP expression by EC subtype and disease group.

Recommended positioning:

- This is not a separate dataset project.
- It should be included as an advanced workflow within the main GSE171524 lung endothelial scRNA-seq repository.
- It demonstrates cell-type-specific integration, RPCA integration, endothelial subtype annotation, and disease-group visualization.
