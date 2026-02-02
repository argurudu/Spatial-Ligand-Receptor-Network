# Spatial Ligand-Receptor Network

The proposed pipeline carries out the following steps of analysis on a given set of annotated spatial RNA-seq samples:
- Segment histological regions for downstream analysis
- Run SCEVAN pipeline to infer copy-number alterations and relative tumor-normal ratios within each subset
- Compute eigengenes for cell type markers within each spatial subset
- Calculate ligand and receptor correlations with cell-type eigengenes, and with direct spatial spots
- Aggregate correlations across multiple samples and subsets

## Instructions for Gathering Subset Coordinates
Spatial coordinates separating histological regions (ex: leading edge, cellular tumor, and infiltrating_tumor) were gathered by hovering over the border between the subsets on an interactive Spatial DimPlot, using the code below:
```R
dim_plot = SpatialDimPlot(seurat_object, interactive = TRUE)
```

After segmenting the main Seurat object into histological regions, the subsets were plotted individually to confirm that spatial spots were correctly divided across their border on the tumor slice.
```R
leading_edge_dim_plot = SpatialDimPlot(leading_edge, interactive = TRUE)
cellular_tumor_dim_plot = SpatialDimPlot(cellular_tumor, interactive = TRUE)
```

All (x, y) coordinates separating regions were recorded for the six samples tested, located under Coordinates-For-Histological-Subsets. 
