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
library(SCEVAN)

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

DATA_PATH = "spatial_samples/<sample>"

run_spatial_pipeline = function(sample_id, data_path, coordinates_files, annotated_image_file, remaining_region, sct_assay = "Spatial", cores = 4) {

  #Load Seurat Object with annotated image
  print("Loading Spatial Data for...", sample_id, "...\n")
  image = Read10X_Image(image.dir = file.path(data_path, "spatial"),
                       image.name = annotated_image_file)
  seurat_object = Load10X_Spatial(data.dir = data_path, image = image)

  #Loop over files with coordinates for subsets
  region_subsets = list()
  for (i in seq_along(coordinates_files)) {
    region_name = names(coordinates_files)[i]

}


#Load Seurat Object with Subset-Annotated Image
