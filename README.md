# Spatial-Ligand-Receptor-Network

The proposed pipeline carries out the following steps of analysis on a given set of annotated spatial RNA-seq samples:
- Segment histological regions for downstream analysis
- Run SCEVAN pipeline to infer copy-number alterations and relative tumor-normal ratios within each subset
- Compute eigengenes for cell type markers within each spatial subset
- Calculate ligand and receptor correlations with cell-type eigengenes, and with direct spatial spots
- Aggregate correlations across multiple samples and subsets

## Instructions for Gathering Subset Coordinates
Spatial coordinates separating histological regions (ex: leading edge, cellular tumor, and infiltrating_tumor) were gathered by hovering over the border between regions on an interaction Spatial DimPlot, using the code below:
```R
dim_plot = SpatialDimPlot(seurat_object, interactive = TRUE)
```
Both x and y coordinates separating regions were recorded for all 6 samples tested, located under Coordinates-For-Histological-Subsets. 
