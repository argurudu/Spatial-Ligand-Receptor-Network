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
library(readr)
install_github("miccec/yaGST")
install_github("AntonioDeFalco/SCEVAN")
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

#' Add SCEVAN metadata to Seurat object
#'
#' Merges selected columns from SCEVAN results into Seurat metadata.
#'
#' @param seurat_obj Seurat object
#' @param results_obj Data frame returned by SCEVAN
#' @param prefix Prefix for metadata column names
#' @return Seurat object with SCEVAN metadata added
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

#' Compute eigengene from expression matrix
#'
#' Performs PCA and returns the first principal component as the eigengene.
#'
#' @param m Expression matrix
#' @return List containing PC1 vector and variance explained
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

#' Compute eigengene vector for a gene set
#'
#' Calculates eigengene values for selected genes.
#'
#' @param seurat_obj Seurat object
#' @param genes Vector of gene names
#' @return Eigengene vector
compute_eigengene_vector = function(seurat_obj, genes) {
  expr_mat = seurat_obj@assays$SCT@scale.data[genes, , drop=FALSE]
  eig_res = get_eigengene(expr_mat)
  return(eig_res[[1]])
}

#' Calculate Spearman correlations
#'
#' Computes Spearman correlations between gene expression and eigengenes
#'
#' @param expr_mat Expression matrix
#' @param eigengenes Eigengene matrix
#' @return Correlation matrix
correlation_calc = function(expr_mat, eigengenes) {
  valid_cols = colSums(is.na(eigengenes)) < nrow(eigengenes)
  cor_valid = cor(t(expr_mat), eigengenes[, valid_cols], method = "spearman")
  cor_full = matrix(NA, nrow = nrow(cor_valid), ncol = ncol(eigengenes),
                    dimnames = list(rownames(cor_valid), colnames(eigengenes)))
  cor_full[, valid_cols] = cor_valid
  return(cor_full)
}

#' Run spatial transcriptomics analysis pipeline
#'
#' Loads Visium data, subsets annotated regions, runs SCEVAN, computes eigengenes,
#' and correlates ligand/receptor expression with cell-type signatures.
#'
#' @param sample_id Sample name
#' @param data_path Path to spatial data
#' @param coordinates_files CSV files defining histological region borders within spatial sample
#' @param annotated_image_file Manually annotated tissue image
#' @param remaining_subset Name for unannotated cells (remaining group not included in coordinates_files)
#' @param sct_assay Assay used for SCTransform
#' @param cores Number of CPU cores
#' @return Writes correlation CSV files
run_spatial_pipeline = function(sample_id, data_path, coordinates_files, annotated_image_file, remaining_subset, eigengene_markers, lr_pairs, sct_assay = "Spatial", cores = 4) {

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
    region_subsets[[region]] = SCTransform(region_subsets[[region]], assay = sct_assay, verbose = FALSE)
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
  ligands = unique(lr_pairs$Ligand)
  receptors = unique(lr_pairs$Receptor)

  region_correlations = list()
  for (region in names(region_subsets)) {
    cat("Computing correlations for region:", region, "...\n")
    expr_mat = region_subsets[[region]]@assays$SCT@scale.data
    eigengenes = region_eigengenes[[region]]
    
    #Ligand correlations
    ligands_present = intersect(ligands, rownames(expr_mat))
    cor_ligands = correlation_calc(expr_mat[ligands_present, , drop=FALSE], eigengenes)
    
    #Receptor correlations
    receptors_present = intersect(receptors, rownames(expr_mat))
    cor_receptors = correlation_calc(expr_mat[receptors_present, , drop=FALSE], eigengenes)
    
    #Store results
    region_correlations[[region]] = list(
      ligand_vs_eigengene = cor_ligands,
      receptor_vs_eigengene = cor_receptors
    )
    write.csv(cor_ligands, paste0("LR_correlations/", sample_id, "_", region, "_ligand_vs_eigengene.csv"))
    write.csv(cor_receptors, paste0("LR_correlations/", sample_id, "_", region, "_receptor_vs_eigengene.csv"))
  }
  return(region_subsets)
}

#' Compute median eigengene correlations across regions
#'
#' Aggregates ligand and receptor correlations across samples and
#' identifies meaningful LR pairs based on correlation thresholds.
#'
#' @param correlation_threshold Minimum correlation cutoff
#' @param histological_regions List of region names
#' @param ligand_dir Directory of ligand correlation files
#' @param receptor_dir Directory of receptor correlation files
#' @param lr_pairs_file Ligand-receptor pair CSV
#' @param output_dir Output directory
#' @return Writes median correlations and filtered LR pairs
median_eigengene_correlations = function(
  correlation_threshold,
  histological_regions,
  ligand_dir,
  receptor_dir,
  lr_pairs,
  output_dir
) {
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  #Read LR pairs
  lr_pairs = lr_pairs %>% distinct(Ligand, Receptor, .keep_all = TRUE)
  
  #Define markers
  non_tumor_eigengenes = c('macrophage_m1','macrophage_m2','bcell','endothelial','smooth_muscle',
                            'cd4_tcell','nk','cd8_tcell','neutrophil','mg_activated','oligo','monocyte',
                            'mg_quiescent','tumor')
  tumor_eigengenes = c('tumor')
  
  #Function to compute median correlation per region
  compute_median = function(files_dir, prefix) {
    region_medians = list()
    for(region in histological_regions) {
      files = list.files(files_dir, pattern = region, full.names = TRUE)
      correlation_list = lapply(files, function(f) read_csv(f, col_types = cols(), show_col_types = FALSE))
      
      #Get all unique rownames and colnames
      all_rows = unique(unlist(lapply(correlation_list, function(df) df[[1]])))
      all_cols = unique(unlist(lapply(correlation_list, function(df) colnames(df)[-1])))
      
      #Build median correlation matrix
      median_mat = matrix(NA, nrow = length(all_rows), ncol = length(all_cols),
                           dimnames = list(all_rows, all_cols))
      for (r in all_rows) {
        for (c in all_cols) {
          vals = sapply(correlation_list, function(df) {
            if (r %in% df[[1]] && c %in% colnames(df)) {
              df_val = df[df[[1]] == r, c, drop = TRUE]
              if (length(df_val) == 0) return(NA) else return(df_val)
            } else return(NA)
          })
          vals = vals[!is.na(vals)]
          if(length(vals) > 0) median_mat[r, c] = median(vals)
        }
      }
      
      region_medians[[region]] = as.data.frame(median_mat)
      write_csv(region_medians[[region]], file.path(output_dir, paste0(prefix, "_vs_eigengene_correlation_spearman_median_", region, ".csv")))
    }
    return(region_medians)
  }
  
  #Compute medians for ligands and receptors
  ligand_medians = compute_median(ligand_dir, "ligand")
  receptor_medians = compute_median(receptor_dir, "receptor")
  
  #Filter meaningful LR pairs
  results = data.frame()
  for(ligand_region in histological_regions) {
    for(receptor_region in histological_regions) {
      ligand_df = ligand_medians[[ligand_region]]
      receptor_df = receptor_medians[[receptor_region]]
      
      filtered_ligands = expand.grid(Ligand = rownames(ligand_df), Sending_Cell_Type = non_tumor_eigengenes) %>%
        rowwise() %>%
        mutate(Ligand_Correlation = ligand_df[ Ligand, Sending_Cell_Type ]) %>%
        filter(!is.na(Ligand_Correlation) & Ligand_Correlation >= correlation_threshold)
      
      filtered_receptors = expand.grid(Receptor = rownames(receptor_df), Receiving_Cell_Type = tumor_eigengenes) %>%
        rowwise() %>%
        mutate(Receptor_Correlation = receptor_df[ Receptor, Receiving_Cell_Type ]) %>%
        filter(!is.na(Receptor_Correlation) & Receptor_Correlation >= correlation_threshold)
      
      meaningful_pairs = lr_pairs %>%
        inner_join(filtered_ligands, by = c("Ligand")) %>%
        inner_join(filtered_receptors, by = c("Receptor")) %>%
        mutate(Ligand_Histological_Type = ligand_region,
               Receptor_Histological_Type = receptor_region) %>%
        select(Ligand_Histological_Type, Ligand, Sending_Cell_Type, Ligand_Correlation,
               Receptor_Histological_Type, Receptor, Receiving_Cell_Type, Receptor_Correlation)
      
      results = bind_rows(results, meaningful_pairs)
    }
  }
  write_csv(results, file.path(output_dir, "meaningful_pairs.csv"))
}

#' Compute direct ligand-receptor correlations for a histological region
#'
#' Calculates Spearman correlations for ligand-receptor pairs within one Seurat subset.
#'
#' @param subset_obj Seurat object subset
#' @param subset_label Name of region/subset
#' @param lr_pairs Ligand-receptor dataframe
#' @return Dataframe of correlations and p-values
compute_direct_lr_correlations = function(subset_obj, subset_label, lr_pairs) {
  
  #Filter LR pairs present in the subset
  lr_pairs_present = lr_pairs %>%
    filter(Ligand %in% rownames(subset_obj) & Receptor %in% rownames(subset_obj))
  if(nrow(lr_pairs_present) == 0) return(NULL)
  
  #Get expression matrix
  expr_mat = as.matrix(GetAssayData(subset_obj, assay = "SCT", slot = "data"))
  
  #Rank transform rows
  expr_rank = t(apply(expr_mat, 1, rank, ties.method = "average"))
  n = ncol(expr_mat)
  
  #Compute correlations for each LR pair
  results = lr_pairs_present %>%
    rowwise() %>%
    mutate(
      rho = cor(expr_rank[Ligand, ], expr_rank[Receptor, ]),
      tval = rho * sqrt((n - 2) / (1 - rho^2)),
      pval = 2 * pt(-abs(tval), df = n - 2),
      subset = subset_label
    ) %>%
    ungroup() %>%
    select(subset, Ligand, Receptor, rho, pval)
  return(results)
}

#' Compute LR correlations across multiple subsets
#'
#' Runs ligand-receptor correlation analysis across multiple regions and
#' summarizes median correlations by histological type.
#'
#' @param subsets_list List of Seurat objects
#' @param subset_names List of subset names
#' @param lr_pairs Ligand-receptor dataframe
#' @return Writes summary CSV file
compute_lr_correlations_multi = function(subsets_list, subset_names, lr_pairs) {
  
  #Compute correlations for all subsets
  all_results = mapply(
    compute_direct_lr_correlations,
    subset_obj = subsets_list,
    subset_label = subset_names,
    MoreArgs = list(lr_pairs = lr_pairs),
    SIMPLIFY = FALSE
  ) %>%
    bind_rows()
  
  #Pivot to wide format and compute median correlations per histological type
  df_wide = all_results %>%
    select(subset, Ligand, Receptor, rho) %>%
    pivot_wider(
      id_cols = c(Ligand, Receptor),
      names_from = subset,
      values_from = rho,
      values_fn = list(rho = mean)
    )
  
  #Compute median correlations for each histological type
  df_wide = df_wide %>%
    rowwise() %>%
    mutate(
      median_leading_edge = median(c_across(starts_with("leading_edge")), na.rm = TRUE),
      median_cellular_tumor = median(c_across(starts_with("cellular_tumor")), na.rm = TRUE),
      median_infiltrating_tumor = median(c_across(starts_with("infiltrating_tumor")), na.rm = TRUE)
    ) %>%
    ungroup()
  
  #Save to file
  write.csv(df_wide, "direct_lr_correlations.csv", row.names = FALSE)
}
