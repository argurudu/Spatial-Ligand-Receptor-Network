library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(readxl)
library(biomaRt)
library(broom)
library(devtools)
install_github("miccec/yaGST")
install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
install_github('plaisier-lab/scSignalMap/scSignalMap')
library(scSignalMap)

#Canonical Markers for all clustered cell types, appended to differentially-expressed cell type markers produced by scSignalMap
MANUALLY_DEFINED_MARKERS = list(
  oligo = c('MBP', 'PLP1','MOG','MOBP','GPR37'),
  mg_activated = c('SPP1','TNF','AIF1'),
  neutrophil = c('CD16','CD15','MPO','CSF3R','ELANE'),
  macrophage_m1 = c('NOS2','IL1B','CXCL10','CD80'),
  macrophage_m2 = c('MRC1','IL10','ARG1','CCL18'),
  monocyte = c('S100A8','S100A9','LYZ','CD163','FCN1'),
  mg_quiescent = c('TMEM119','P2RY12','CX3CR1'),
  nk = c('NCAM1','KLRG1','NCR1','IL15RA','IFNG'),
  bcell = c('CD19','CD22','CD40','CD38','CXCR5'),
  endothelial = c('PECAM1','CLDN5','ICAM2','CDH5','ESAM','VWF'),
  smooth_muscle = c('ACTA2','CNN1','CALD1','MYH11','TAGLN','TAGLN2'),
  cd4_tcell = c('CD4','CD3D','CD40LG'),
  cd8_tcell = c('CD8A','CD8B','CD27','CD28','PRF1','PDCD1'),
  tumor = c('PRC1','CDKN2A', 'PTPRZ1')
)

add_scevan_meta = function(seurat_obj, results_obj, prefix = "SCEVAN") {
  if (is.data.frame(results_obj) && !is.null(rownames(results_obj))) {
    common = intersect(rownames(results_obj), colnames(seurat_obj))
    if (length(common) == 0) return(seurat_obj)
    
    want_cols = c("subclone", "clone", "SCEVAN_clone", "malignancy", "aneuploidy_score", "status")
    present = intersect(want_cols, colnames(results_obj))
    if (length(present) == 0) present = colnames(results_obj)
    
    for (col in present) {
      vect = rep(NA, ncol(seurat_obj))
      names(vect) = colnames(seurat_obj)
      vect[common] = as.character(results_obj[common, col])
      seurat_obj[[paste0(prefix, "_", col)]] = vect
    }
    return(seurat_obj)
  }
  return(seurat_obj)
}

get_eigengene = function(m) {
  pca1 = prcomp(m)
  varExp = summary(pca1)$importance[2,1]
  pc1 = pca1$rotation[,'PC1']
  if(sum(cor(pc1,t(m))) > 0) {
    return(list(pc1, varExp))
  } else {
    return(list(-pc1,varExp))
  }
}

compute_eigengene_vector = function(seurat_obj, genes) {
  expr_mat = seurat_obj@assays$SCT@scale.data[genes, , drop=FALSE]
  eig_res = get_eigengene(expr_mat)
  return(eig_res[[1]])
}

correlation_calc = function(expr_mat, eigengenes) {
  valid_cols = colSums(is.na(eigengenes)) < nrow(eigengenes)
  cor_valid = cor(t(expr_mat), eigengenes[, valid_cols], method = "spearman")
  cor_full = matrix(NA, nrow = nrow(cor_valid), ncol = ncol(eigengenes),
                    dimnames = list(rownames(cor_valid), colnames(eigengenes)))
  cor_full[, valid_cols] = cor_valid
  return(cor_full)
}

##############################
#### Main Seurat Pipeline ####
##############################

#Run scSignalMap on cell-clustered, scRNA-seq GBM sample to produce differentially-expressed markers and LR pairs
scRNA_object = readRDS("/Users/SRG15/Desktop/spatial-seq/GSE197543/GSE197543_normalized_ensembl.rds")
interactions = MapInteractions(
  seurat_obj,
  group_by = "CellAnnotationMerged",
  avg_log2FC_gte = 0.25,
  p_val_adj_lte = 0.05,
  min_pct = 0.1,
  species = "human",
  gene_id = "ensembl"
)
write.csv(interactions, "GSE197543_scSignalMap.csv")

#Function to run spatial analysis pipeline after scSignalMap
run_spatial_pipeline = function(sample_id, data_path, coordinates_files, annotated_image_file, remaining_subset, sct_assay = "Spatial", cores = 4) {

  #Load Seurat object with annotated image
  cat("Loading Spatial Data for...", sample_id, "...\n")
  image = Read10X_Image(image.dir = file.path(data_path, "spatial"),
                       image.name = annotated_image_file)
  seurat_object = Load10X_Spatial(data.dir = data_path, image = image)

  #Loop over files with coordinates for subsets
  region_subsets = list()
  for (i in seq_along(coordinates_files)) {
    filename = basename(coordinates_files[i])
    region_name = sub("_[^_]+_coordinates\\.csv$", "", filename)
    coordinates_file = coordinates_files[[i]]
    coordinates_df = read.csv(coordinates_file)

    #Create segmentation
    segmentation = CreateSegmentation(coordinates_df)
    seurat_object[[region_name]] = Overlay(seurat_object[["slice1"]], segmentation)

    #Subset the Seurat object
    region_subsets[[region_name]] = subset(seurat_object, cells = Cells(seurat_object[[region_name]]))
  }

  #Create the remaining subset
  all_annotated_cells = unlist(lapply(region_subsets, Cells))
  region_subsets[[remaining_subset]] = subset(seurat_object, cells = setdiff(Cells(seurat_object), all_annotated_cells))

  #Run SCEVAN on all subsets
  scevan_results = list()
  for (region in names(region_subsets)) {
    cat("Running SCEVAN on region:", region, "...\n")
    counts = as.matrix(GetAssayData(region_subsets[[region]], assay = "Spatial", slot = "counts"))
    results = pipelineCNA(counts,
                      sample = paste0(sample_id, "_", region),
                      par_cores = cores,
                      SUBCLONES = FALSE,
                      beta_vega = 0.5)
    scevan_results[[region]] = results
  }

  #Add SCEVAN results to metadata
  for (region in names(region_subsets)) {
    region_subsets[[region]] = add_scevan_meta(region_subsets[[region]], scevan_results[[region]])
  }

  #Run SCTransform
  for (region in names(region_subsets)) {
    cat("Running SCTransform on region: ", region, "...")
    region_subsets[[region]] <- SCTransform(region_subsets[[region]], assay = sct_assay, verbose = FALSE)
  }

  #Compute eigengenes for each subset
  region_eigengenes = list()
  for (region in names(region_subsets)) {
    cat("Computing eigengenes for region:", region, "...\n")
    region_eigengenes[[region]] = sapply(names(eigengene_markers), function(ct) {
      all_genes = eigengene_markers[[ct]]
      present_genes = all_genes[all_genes %in% rownames(region_subsets[[region]]@assays$SCT@scale.data)]
      frac = length(present_genes) / length(all_genes)
      if (frac < 0.5) {
        return(rep(NA, ncol(region_subsets[[region]]@assays$SCT@scale.data)))
      }
      compute_eigengene_vector(region_subsets[[region]], present_genes)
    }, simplify = FALSE)
    
    region_eigengenes[[region]] = do.call(cbind, region_eigengenes[[region]])
    colnames(region_eigengenes[[region]]) = names(eigengene_markers)
  }

  #Compute correlations between ligands/receptors and eigengenes
  lr_pairs = read.csv('/Users/SRG15/Desktop/spatial-seq/LR/LR_pairs_output_GSE197543.csv', header=TRUE)
  ligands = unique(lr_pairs$Ligand)
  receptors = unique(lr_pairs$Receptor)

  region_correlations = list()
  for (region in names(region_subsets)) {
    cat("Computing correlations for region:", region, "...\n")
    expr_mat = region_subsets[[region]]@assays$SCT@scale.data
    eigengenes = region_eigengenes[[region]]
    
    # Ligand correlations
    ligands_present = intersect(ligands, rownames(expr_mat))
    cor_ligands = correlation_calc(expr_mat[ligands_present, , drop=FALSE], eigengenes)
    
    # Receptor correlations
    receptors_present = intersect(receptors, rownames(expr_mat))
    cor_receptors = correlation_calc(expr_mat[receptors_present, , drop=FALSE], eigengenes)
    
    # Store results
    region_correlations[[region]] = list(
      ligand_vs_eigengene = cor_ligands,
      receptor_vs_eigengene = cor_receptors
    )
  }
  write.csv(cor_ligands, paste0("LR_correlations/", sample_id, "_", region, "_ligand_vs_eigengene.csv"))
  write.csv(cor_receptors, paste0("LR_correlations/", sample_id, "_", region, "_receptor_vs_eigengene.csv"))
}

#Calculate median eigengene correlations
median_eigengene_correlations = function(){
  
}
