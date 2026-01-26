# Spatial-Ligand-Receptor-Network

The proposed pipeline carries out the following steps of analysis on a given set of annotated spatial RNA-seq samples:
- Segment histological regions for downstream analysis
- Run SCEVAN pipeline to infer copy-number alterations and relative tumor-normal ratios within each subset
- Compute eigengenes for cell type markers within each spatial subset
- Calculate ligand and receptor correlations with cell-type eigengenes, and with direct spatial spots
- Aggregate correlations across multiple samples and subsets

```R
print('hello world')
```
