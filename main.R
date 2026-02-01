##############################
#### Main Seurat Pipeline ####
##############################

#Import functions
source("LigandReceptorNetwork.R")

#Run scSignalMap on cell-clustered scRNA-seq GBM sample to produce differentially-expressed markers and LR pairs
scRNA_object = readRDS("scRNA-seq_sample.rds")
interactions = MapInteractions(
  scRNA_object,
  group_by = "CellAnnotationMerged", #name of column with cell type annotations
  avg_log2FC_gte = 0.25,
  p_val_adj_lte = 0.05,
  min_pct = 0.1,
  species = "human",
  gene_id = "ensembl"
)
write.csv(interactions, "scSignalMap_results.csv")

#spatialRNA-seq Sample Names
samples = c(
  "UKF_251",
  "UKF_243",
  "UKF_260",
  "UKF_266",
  "UKF_269",
  "UKF_334"
)

for (sample_id in samples) {
  cat("Processing sample:", sample_id, "\n")
  coord_file = paste0("leading_edge_", sample_id, "_coordinates.csv")
  run_spatial_pipeline(
    sample_id,
    paste0(sample_id, "_data/spatial"),
    c(coord_file),
    paste0(sample_id, "/spatial/tissue_lowres_image_annotated.png"),
    "cellular_tumor",
    sct_assay = "Spatial",
    cores = 4
  )
}
#Compute medians after all samples run
median_eigengene_correlations(
  0.3,
  histological_regions = c("leading_edge", "cellular_tumor", "infiltrating_tumor"),
  ligand_dir = "L_vs_markers",
  receptor_dir = "R_vs_markers",
  lr_pairs_file = "scSignalMap_results.csv",
  output_dir = "mined_correlations"
)

#Compute direct ligand-receptor correlations on all samples
lr_pairs = read.csv("scSignalMap_results.csv")
compute_lr_correlations_multi(
  subsets_list,
  subset_names,
  lr_pairs
)
