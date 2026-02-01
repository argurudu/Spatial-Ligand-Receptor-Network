##############################
#### Main Seurat Pipeline ####
##############################

#Import functions
source("LigandReceptorNetwork.R")

#Run scSignalMap on cell-clustered scRNA-seq GBM sample to produce differentially-expressed markers and LR pairs
scRNA_object = readRDS("/Users/SRG15/Desktop/spatial-seq/GSE197543/GSE197543_normalized_ensembl.rds")
interactions = MapInteractions(
  seurat_obj,
  group_by = "CellAnnotationMerged", #name of column with cell type annotations
  avg_log2FC_gte = 0.25,
  p_val_adj_lte = 0.05,
  min_pct = 0.1,
  species = "human",
  gene_id = "ensembl"
)
write.csv(interactions, "scSignalMap_results.csv")

