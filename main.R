##############################
#### Main Seurat Pipeline ####
##############################

#Import statements
source("LigandReceptorNetwork.R")
install_github('plaisier-lab/scSignalMap/scSignalMap')
library(scSignalMap)

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
write.csv(interactions, "scSignalMap_LR_pairs.csv")

#Loading sample coordinates
sample_coords = list(
  UKF_251 = c(
    "leading_edge_UKF_251_coordinates.csv",
    "infiltrating_tumor_UKF_251_coordinates.csv"
  ),
  UKF_243 = c("Coordinates-For-Histological-Subsets/leading_edge_UKF_243_coordinates.csv"),
  UKF_260 = c("Coordinates-For-Histological-Subsets/leading_edge_UKF_260_coordinates.csv"),
  UKF_266 = c("Coordinates-For-Histological-Subsets/leading_edge_UKF_266_coordinates.csv"),
  UKF_269 = c("Coordinates-For-Histological-Subsets/leading_edge_UKF_269_coordinates.csv"),
  UKF_334 = c("Coordinates-For-Histological-Subsets/leading_edge_UKF_334_coordinates.csv")
)

eigengene_markers = read.csv("eigengene_markers.csv") #file with MANUALLY_DEFINED_MARKERS appended to differentially-expressed markers produced by scSignalMap
lr_pairs = read.csv("scSignalMap_LR_pairs.csv") #file with LR pairs produced by scSignalMap
subsets_list = list()

for (sample_id in names(sample_coords)) {
  cat("Processing sample:", sample_id, "\n")
  region_subsets = run_spatial_pipeline(
    sample_id,
    paste0(sample_id, "_data/spatial"),
    sample_coords[[sample_id]],
    paste0(sample_id, "/spatial/tissue_lowres_image_annotated.png"),
    "cellular_tumor",
    eigengene_markers,
    lr_pairs,
    sct_assay = "Spatial",
    cores = 4
  )
  for (region in names(region_subsets)) {
    subsets_list[[paste0(sample_id, "_", region)]] = region_subsets[[region]]
  }
}
#Compute medians after all samples run
median_eigengene_correlations(
  0.3,
  histological_regions = c("leading_edge", "cellular_tumor", "infiltrating_tumor"),
  ligand_dir = "L_vs_markers",
  receptor_dir = "R_vs_markers",
  lr_pairs_file = "scSignalMap_LR_pairs.csv",
  output_dir = "mined_correlations"
)

#Compute direct ligand-receptor correlations on all samples
lr_pairs = read.csv("scSignalMap_results.csv")
subset_names = c(
  "UKF_251_leading_edge","UKF_251_cellular_tumor","UKF_251_infiltrating_tumor",
  "UKF_243_leading_edge","UKF_243_cellular_tumor"
  "UKF_260_leading_edge","UKF_260_cellular_tumor"
  "UKF_266_leading_edge","UKF_266_cellular_tumor"
  "UKF_269_leading_edge","UKF_269_cellular_tumor"
  "UKF_334_leading_edge","UKF_334_cellular_tumor"
)
compute_lr_correlations_multi(
  subsets_list,
  subset_names,
  lr_pairs
)
