
## ---------------------------
##
## Script name: seurat_utils.R
##
## Purpose of script: # Helper functions to work with Seurat objects 
##
## Author: Ann Strange
##
## Date Created: 2024-04-03
##
## Email: ann.strange@cuanschutz.edu
##
## ---------------------------
##
## Notes:
## Unless otherwise specified these functions are designed to work with Seurat v5 objects   
##
## ---------------------------

## load up our functions into memory

# source("R/seurat_utils.R")") 

## ---------------------------

require(tidyverse)
require(data.table)
library(readxl)
library(Seurat)

options(Seurat.object.assay.version = "v5")

# Get metadata from Seurat object
get_metadata <- function(seurat_obj) {
  metadata <- seurat_obj@meta_data 
  metadata <- as.data.frame(metadata)
  return(metadata)
}

#' Load a Seurat object from a RDS file
#'
#' @param seurat.filename The filename of the Seurat .RDS object
#' @return A Seurat object
#' @examples
#' sobj <- load_seurat_object("data/seuratObject.RDS")
#' sobj <- load_seurat_object(paste0(psr.data.dir, "d2535822-e266-4373-817a-6bb17917b354_seuratObject.RDS")
load_seurat_object <- function(seurat.filename) {
  #psr.seurat.file <- paste0(psr.data.dir, "d2535822-e266-4373-817a-6bb17917b354_seuratObject.RDS")
  sobj <- readRDS(seurat.filename)
  Assays(sobj)
  # print some basic info
  print(paste0("Number of cells: ", nrow(sobj)))
  print(paste0("Number of genes: ", ncol(sobj)))
  print(paste0("Number of assays: ", length(Assays(sobj))))
  # print the metadata column names
  print(colnames(sobj@meta.data))
  print(Layers(sobj)) # counts data
  print(Assay(sobj)) # RNA is a common default
  return(sobj)
}


#' With a Seurat object created with LoadNanostring(), add the metatdata from a flat file
#' 
#' @param sobj The Seurat object
#' @param metadata_file The metadata file to load
#' @return The Seurat object with metadata rownames relabeled
#' @examples
#' sobj <- rename_keys(sobj)
load_mouse_meta_load_missed <- function(sobj, metadata_file) {
  # Seurat LoadNanostring() workaround to load metadata flat file including making sure the keys match
  # this loads morph markers, measurements, X, Y, (global and local) and cell_ID
  
  df_obj_meta <- read.csv(metadata_file)
  #print(colnames(df_obj_meta))
  
  df_obj_meta$cell_num = df_obj_meta$cell_ID 
  df_obj_meta$key = paste0( df_obj_meta$cell_ID, "_", df_obj_meta$fov)
  rownames(df_obj_meta) <- df_obj_meta$key
  
  head(rownames(sobj@meta.data))
  print("Number of cells in Seurat object before adding metadata")
  print(nrow(sobj@meta.data))
  # AddMetaData needs the same records in the same order
  print("Do the number of rows in the metadata file match the sobj?")
  print(nrow(df_obj_meta) == nrow(sobj@meta.data)) # FALSE
  lost_on_load <- setdiff(rownames(df_obj_meta), rownames(sobj@meta.data)) # about 100, drop these 
  length(lost_on_load)
  if(length(lost_on_load) > 0){
    print ("These are the rows that are in the metadata file but not the Seurat obj, we will drop them from the metadata file.")
    print(head(lost_on_load))
  }
  #print("expect 0 cell differences")
  #print(setdiff(rownames(sobj@meta.data), rownames(df_obj_meta))) # 0
  
  # delete these from df_obj_metadata
  df_obj_meta <- df_obj_meta[!rownames(df_obj_meta) %in% lost_on_load,]
  #print(setdiff(rownames(df_obj_meta), rownames(sobj@meta.data))) # should now be 0
  
  df_obj_meta$cell_ID <- df_obj_meta$key # restore cell_ID in standard format
  sobj <- AddMetaData(sobj, metadata = df_obj_meta)
  #head(sobj@meta.data)
  #head(rownames(sobj@meta.data))
  
  return (sobj)  
}


#' With a Seurat object, rename cells to be 'cell_id' in the metadata
#' # This rename is needed as LoadNanostring() uses <cell>_<fov> cell IDs vs the c_<slide>_<fov>_<cell> 
#' format used elsewhere in the cosMx pipeline
#' 
#' @param sobj The Seurat object
#' @return The Seurat object with metadata rownames relabeled
#' @examples
#' sobj <- rename_keys(sobj)
rename_keys <- function(sobj) {
  # given a Seurat object, rename cells to be 'cell_id' in the metadata
  
  print("Renaming keys for cells in Seurat object to match cosMx format")
  new_cells <- sobj@meta.data$cell
  # check format
  print(paste("Example of new cell rownames", paste(new_cells[1:5], collapse=' ')))
  sobj <- RenameCells(sobj, new.names = new_cells)
  
  # spot check rename applied in another place:
  #print(paste("Example of new cell rownames in colnames", paste(head(colnames(sobj),3), collapse=' ')))

  return (sobj)
}


#' Compare cell IDs in Seurat object to the dataframe
#' This is useful to check before calling AddMetaData to make sure all the keys match up.
#'
#' @param sobj The Seurat object
#' @param df The dataframe for the metadata file
#' @return success boolean, True if AddMetaData can be called successfully
#' @examples
#' success <- do_rownames_match(seurat_obj, df_metadata)
do_rownames_match <- function(sobj, df) {

  keys1 <- rownames(sobj@meta.data)
  keys2 <- rownames(df)
  if (length(keys1) != length(keys2)) {
    print("Number of cells don't match. Check the keys.")
  }
  print(paste("Number of cells in Seurat object before adding metadata", length(keys1)))
  print(paste("Number of cells in metadata file", length(keys2)))
  
  diff1 <- setdiff(rownames(keys1), rownames(keys2))
  diff2 <- setdiff(rownames(keys2), rownames(keys1))
  if (length(diff1) > 0) {
    print("These keys are in the Seurat object but not the metadata file")
    print(diff1)
  }
  if (length(diff2) > 0) {
    print("These keys are in the metadata file but not the Seurat object")
    print(diff2)
  }
  if (length(diff1) == 0 & length(diff2) == 0) {
    print("All keys match")
    if (identical(keys1, keys2)) {
      print("The keys are in the same order")
      return (TRUE)
    } else {
      print("The keys match, but reorder.")
      return (FALSE)
    }
  }
  return (FALSE)
}
  

sanitize_name <- function(marker){
  # remove special characters from marker name; replace with -
  # good for marker that might contain "/" and filename text that might contain special characters
  return(gsub("[^[:alnum:]]", "_", marker))
}
#sanitize_marker_name("CCL3/L1")

# get percent positive for a marker
get_percent_pos <- function(df_meta, marker){
  formatted_pct <- sprintf("%.2f%%", sum(df_meta[[marker]]) / length(df_meta[[marker]]) * 100)
  return(formatted_pct)
}
#get_percent_pos(df_metadata, "CD3E")



#' subtract the negative controls from the counts matrix
#' and return the modified sparse count matrix, 0 values in the Negative features 
#' Assumes the Negative features are still present in the Seurat object 
#' and the Assay is "Nanostring", layer is "counts", per defaults from Seurat::LoadNanostring()
#'
#' @param obj The Seurat object
#' @return modified count matrix
#' @examples
#' obj <- subtract_neg_control(obj)
subtract_neg_control <- function(obj, assay="RNA", layer="data"){
  # get negative probe feature indexes, pre-normalization 
  counts_matrix <- GetAssayData(object = obj, assay = assay, layer=layer) ## or should use normalized?
  
  feature_names <- rownames(obj)
  # which feature_names in list begin with "Neg"?
  
  feature_names 
  # Filter feature names that begin with "Neg"
  neg_feature_names <- feature_names[grep("^Negative", feature_names)]
  #print(neg_feature_names)
  
  # Get row indices for neg_feature_names
  neg_indices <- which(feature_names %in% neg_feature_names)
  non_neg_indices <- which(!(feature_names %in% neg_feature_names))
  # Create sparse matrices for neg_feature_names and others
  neg_matrix <- counts_matrix[neg_indices, ]
  non_neg_matrix <- counts_matrix[non_neg_indices, ]
  
  # Verify the results
  print(dim(neg_matrix))       # Dimensions of the matrix with neg_feature_names
  #print(length(sys_control_names))
  print(dim(non_neg_matrix))   # Dimensions of the matrix without neg_feature_names
  
  filler_matrix <- Matrix(0, nrow(neg_matrix), ncol(neg_matrix),  sparse = TRUE)
  rownames(filler_matrix) <- c(neg_feature_names)
  
  negmean <- neg_matrix %>%
    as.matrix() %>%
    t() %>%
    Matrix::rowMeans()
  head(negmean)
  
  # subtract cell-wise negative 'background' from counts matrix
  count_replace <- sweep(non_neg_matrix, 2, negmean, FUN = "-")
  
  #count_replace <- apply(non_neg_matrix,1,function(row) row - negmean) # mem hog, don't! 
  count_replace[count_replace < 0 ] <- 0
  count_replace <- ceiling(count_replace)
  
  # replace the matrix
  # Assume new_counts is your new matrix and it's already properly formatted
  # Check that the dimensions and names match
  if (all(rownames(count_replace) == rownames(non_neg_matrix )) &&
      all(colnames(count_replace) == colnames(non_neg_matrix))) {
    
  } else {
    stop("Error: Dimension names of new_counts do not match those in seurat_object")
  }
  
  print(paste("Negative control subtraction complete.", sum(non_neg_matrix) - sum(count_replace),
              " counts removed of ", sum(non_neg_matrix), " total counts."))
  print(paste("Total negative probe counts: ", sum(neg_matrix), " so any difference is due to rounding, or neg probes > pos per cell."))
  
  # restore lost features with 0s
  count_replace <- rbind(count_replace, filler_matrix)
  # SetAssayData(obj_list[[i]], layer = "counts", new.data = remove_sys_control(obj_list[[i]], "Nanostring", "counts"))
  obj <- SetAssayData(object = obj, assay = assay, layer=layer, new.data = count_replace)
  
  return(obj)

}

#' remove the system control features from the counts matrix
#' and return the modified Seurat object.  
#'
#' @param obj The Seurat object
#' @param assay The assay to use, default is "Nanostring"
#' @param layer The layer to use, default is "counts"
#' @return sobj with system control features removed
#' @examples
#' obj <- remove_sys_control(obj)
#' obj <- remove_sys_control(obj, assay="Nanostring", layer="counts")
remove_sys_control <- function(obj, assay="RNA", layer="counts"){
  # get negative probe feature indexes, pre-normalization 
  counts_matrix <- GetAssayData(object = obj, assay = assay, layer=layer) 
  
  feature_names <- rownames(obj)
  # which feature_names in list begin with "Neg"?
  
  feature_names 
  sys_control_names <- feature_names[grep("^SystemControl", feature_names)]
  genes_to_keep <- feature_names[!feature_names %in% sys_control_names]
  
  obj <- subset(obj, features = genes_to_keep)
  print(paste("Remaining dimensions:", paste(dim(obj), collapse = " x ")))
  
  return(obj)
  
}

#' retrieve the cell-wise mean of the negative control features
#' Assumes: Negative control features are still present in the Seurat object
#'
#' @param obj The Seurat object
#' @return numeric vector of mean negative control values by cell_id
#' @examples
#' obj <- get_neg_control_means(obj)
get_neg_control_means <- function(obj, assay="RNA", layer="counts"){
  
  counts_matrix <- GetAssayData(object = obj, assay = assay, layer=layer) ## or should use normalized?
  
  feature_names <- rownames(obj)
  # which feature_names in list begin with "Neg"?
  
  feature_names 
  neg_feature_names <- feature_names[grep("^Negative", feature_names)]
  #print(neg_feature_names)
  
  # Get row indices for neg_feature_names
  neg_indices <- which(feature_names %in% neg_feature_names)
  # Create sparse matrices for neg_feature_names and others
  neg_matrix <- counts_matrix[neg_indices, ]
  
  negmean <- neg_matrix %>%
    as.matrix() %>%
    t() %>%
    Matrix::rowMeans()
  head(negmean)
  
  # throw error if negmeans are all zero
  if(all(negmean == 0)){
    stop("Error: All negative control means are zero. Check that negative control features are present in the counts matrix.")
  }
  
  return(negmean)
  
}

#' Given a nanostring/Bruker provided gene annotations file, extract the key bits in a single table.
#'
#' @param gene_annotations_file e.g. LBL-11178-03-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx 
#' @return dataframe of cosMx gene display names with annotations
#' @examples
#' df <- get_nanostring_gene_annotations ("path/to/gene_annotations.xlsx")
get_nanostring_mouse_gene_annotations <- function (gene_annotations_file) {
  
  df_gene_data <- as.data.frame(read_excel(ns_gene_fn, sheet = 2, skip = 1, col_names = TRUE))
  head(df_gene_data)
  print(nrow(df_gene_data))
  colnames(df_gene_data) <- sanitize_name(colnames(df_gene_data))
  # remove copywrite line by Gene_Symbol_s_ not NA
  df_gene_data <- df_gene_data[!is.na(df_gene_data$Gene_Symbol_s_),]
  df_gene_data <- df_gene_data[,c("Display_Name", "Human_Gene", "Gene_Symbol_s_", "Gene_Name_s_", "Alias_es_" )]
  df_gene_data <- df_gene_data %>% dplyr::rename(Display_Name = Display_Name, 
                                                 Gene_Symbol = Gene_Symbol_s_, 
                                                 Gene_Long_Name = Gene_Name_s_, 
                                                 Aliases = Alias_es_)
  
  df_gene_annots <- as.data.frame(read_excel(ns_gene_fn, sheet = 3, skip = 1, col_names = TRUE))
  head(df_gene_annots)
  print(nrow(df_gene_annots))
  df_gene_annots <- df_gene_annots[, !names(df_gene_annots) %in% c("Core", "Add-On")]
  colnames(df_gene_annots) <- sanitize_name(colnames(df_gene_annots))
  df_gene_annots$Annotations <- apply(df_gene_annots[, -1], 1, function(x) {
    # Get the column names where the value is "+"
    paste(names(df_gene_annots)[-1][x == "+"], collapse = ", ")
  })
  head(df_gene_annots$Annotations)
  # merge df_gene_data and df_gene_annots by Gene
  df_gene_annots <- df_gene_annots[,c("Gene", "Annotations")]
  df_gene_data <- merge(df_gene_data, df_gene_annots, by.x = "Display_Name", by.y = "Gene")
  head(df_gene_data)
  
  return(df_gene_data)
  
}

#' Given a nanostring/Bruker provided gene annotations file, extract the list of genes annotated for cell typing.
#'
#' @param gene_annotations_file e.g. LBL-11178-03-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx 
#' @return list of genes with a "+" in the requested annotation column
#' @examples
#' genes_for_celltyping <- get_nanostring_celltyping_genes ("path/to/gene_annotations.xlsx", "Cell Typing")
get_nanostring_celltyping_genes <- function (gene_annotations_file, colname = "Cell Typing") {
  
  df_gene_data <- as.data.frame(readxl::read_excel(ns_gene_fn, sheet = 3, skip = 1, col_names = TRUE))
  head(df_gene_data, 3)
  print(nrow(df_gene_data))
  cols_to_keep <- c("Gene", colname)
  df_gene_data <- df_gene_data %>%
    dplyr::select(all_of(cols_to_keep)) %>%
    dplyr::filter(df_gene_data[[colname]] == "+")
  
  gene_list <- df_gene_data$Gene
  
  return(gene_list) # over 200 genes
  
}


#' Given a dataframe including SampleID and Run_Tissue, return the bounding box of the SampleID
#'
#' @param df_meta dataframe including SampleID, Run_Tissue, and global pixel coordinates
#' @return matrix of x,y min and max coordinates (bounding box)
#' @examples
#' df <- get_bbox_of_sample (sobj@meta.data, "Slide2name", "S30")
get_bbox_of_sample <- function(df_meta, Run_Tissue, SampleID) {
  
  df_meta2 <- df_meta %>%
    dplyr::filter(Sample.ID == SampleID & Run_Tissue_name == Run_Tissue)
  
  padding <- 0
  x_min <- round(min(df_meta2$CenterX_global_px)) - padding
  y_min <- round(min(df_meta2$CenterY_global_px)) - padding
  x_max <- round(max(df_meta2$CenterX_global_px)) + padding
  y_max <- round(max(df_meta2$CenterY_global_px)) + padding
  print(paste(y_min, y_max, x_min, x_max))
  
  #sobj_fov <- subset(sobj, subset = Sample.ID %in% sample_ids )
  # Create a 2x2 matrix
  bounding_box <- matrix(c(y_min, y_max, x_min, x_max), nrow = 2, byrow = TRUE,
                         dimnames = list(c("y", "x"), c("min", "max")))
  
  return(bounding_box)
}

#' Given a dataframe including fov and Run_Tissue, return the bounding box of the fov
#'
#' @param df_meta dataframe including SampleID, Run_Tissue, and global pixel coordinates
#' @param Run_Tissue Run_Tissue_name 
#' @fov_num fov number
#' @return matrix of x,y min and max coordinates (bounding box)
#' @examples
#' df <- get_bbox_of_sample (sobj@meta.data, "Slide2name", "S30")
get_bbox_of_fov <- function(df_meta, Run_Tissue, fov_num) {
  
  df_meta2 <- df_meta %>%
    dplyr::filter(fov == fov_num & Run_Tissue_name == Run_Tissue)
  
  padding <- 0
  x_min <- round(min(df_meta2$CenterX_global_px)) - padding
  y_min <- round(min(df_meta2$CenterY_global_px)) - padding
  x_max <- round(max(df_meta2$CenterX_global_px)) + padding
  y_max <- round(max(df_meta2$CenterY_global_px)) + padding
  print(paste(y_min, y_max, x_min, x_max))
  
  #sobj_fov <- subset(sobj, subset = Sample.ID %in% sample_ids )
  # Create a 2x2 matrix
  bounding_box <- matrix(c(y_min, y_max, x_min, x_max), nrow = 2, byrow = TRUE,
                         dimnames = list(c("y", "x"), c("min", "max")))
  
  return(bounding_box)
}

pixels_to_microns <- function(min_px, max_px) {
  # length in microns
  # an fov is 510 microns^2 and ___ pixels^2.
  factor = 4225 / 510  # px in fov x axis / microns  
  dist_microns = abs(max_px - min_px) / factor
  return (dist_microns)
}



#' Given a Seurat object, return the counts matrix as a dataframe by cell
#'
#' @param sobj Seurat object 
#' @return counts matrix as a dataframe by cell
#' @examples
#' df <- get_counts_df (sobj_subset)
get_counts_df <- function(sobj){
  
  exp_data <- as.matrix(GetAssayData(sobj, assay="RNA", layer="counts"))
  class(exp_data)
  dim(exp_data) # 1000 x 191K
  # how to subset matrix by gene
  
  return(as.data.frame(t(exp_data)))
}



# Summarize the counts matrix by sample
get_sample_pcts <- function(df_counts_mtx, sample_id) {
  
  df_pct_sample <- df_counts_mtx %>%
    dplyr::select(c('Cdkn1a','Mki67','Lmna')) %>%
    summarize(count = n(), # number of cells
              Cdkn1a_count = sum(Cdkn1a > 0),
              Mki67_count = sum(Mki67 > 0),
              Lmna_count = sum(Lmna > 0), 
              .groups = "drop") %>%
    dplyr::filter(count > 0)  %>%
    mutate( sample=sample_id,
            pct_cells_Cdkn1a = (Cdkn1a_count / sum(count)) * 100,
            pct_cells_Mki67 = (Mki67_count / sum(count)) * 100,
            pct_cells_Lmna = (Lmna_count / sum(count)) * 100) %>%
    dplyr::relocate(sample) # move sample first
  
  return(df_pct_sample)
}

# Summarize the counts matrix by sample. Hardcoded colnames, sorry.
get_sample_pcts_bycelltype <- function(df_ct, sample_id, celltype, marker_name) {
  
  df_pct_sample <- df_ct %>%
    dplyr::select(c('Sample.ID', marker_name, 'CellType_Integrated')) %>%
    dplyr::mutate(marker = df_ct[[marker_name]]) %>% 
    dplyr::filter(Sample.ID == sample_id & CellType_Integrated == celltype) %>%
    summarize(cell_count = n(), # number of cells
              pos_for_marker_count = count(marker > 0),
              num_probes = sum(marker),
              .groups = "drop") %>%
    #dplyr::filter(count > 0)  %>%
    mutate(sample=sample_id,
                       marker = marker_name,
                       cell_type = celltype,
                       #num_cells_pos = marker_count,
                       pct_cells_pos = (pos_for_marker_count / cell_count) * 100) %>%
    dplyr::relocate(sample) # move sample first
  
  return(df_pct_sample)
}

# Has some hardcoded columns (sorry)
get_pct_fov_celltype <- function(df_ct, sample_id, cell_type, marker_name){
  
  df_pct_fov_celltype <- df_ct %>%
    dplyr::select(c('fov', all_of(marker_name), 'CellType_Integrated', 'Sample.ID')) %>%
    dplyr::mutate(marker = df_ct[[marker_name]]) %>% 
    dplyr::filter(Sample.ID == sample_id & CellType_Integrated == cell_type) %>%
    group_by(fov) %>%
    summarize(cell_count = n(), # number of cells
              pos_for_marker_count = count(marker > 0),
              num_probes = sum(marker),
              .groups = "drop") %>%
    #dplyr::filter(count > 0)  %>%
    mutate(sample=sample_id,
           marker = marker_name,
           cell_type = cell_type,
           pct_cells_pos = (pos_for_marker_count / cell_count) * 100) %>%
    dplyr::relocate(sample) 
  
  return(df_pct_fov_celltype)
}


# Get Percentage of cells for each fov that are positive, independent of celltype
get_pct_fov <- function(df_ct, sample_id, marker_name){
  
  df_pct_fov <- df_ct %>%
    dplyr::select(c('fov', marker_name, 'CellType_Integrated', 'Sample.ID')) %>%
    dplyr::mutate(marker = df_ct[[marker_name]]) %>% 
    dplyr::filter(Sample.ID == sample_id) %>%
    group_by(fov) %>%
    summarize(cell_count = n(), # number of cells
              pos_for_marker_count = count(marker > 0),
              num_probes = sum(marker),
              .groups = "drop") %>%
    #dplyr::filter(count > 0)  %>%
    mutate(sample=sample_id,
           marker = marker_name,
           pct_cells_pos = (pos_for_marker_count / cell_count) * 100) %>%
    dplyr::relocate(sample) 
  
  return(df_pct_fov)
}


parse_third_token <- function(string) {
  # Split the string by underscores and extract the third token
  tokens <- strsplit(string, "_")[[1]]
  return(tokens[3])
}

# Summarize the counts matrix by fov
get_fov_pcts <- function(df_counts_mtx, sample_id) {
  df_counts_mtx$fov <- sapply(df_counts_mtx$cell, parse_third_token)
  df_pct_fov <- df_counts_mtx %>%
    dplyr::select(c('fov','Cdkn1a','Mki67','Lmna')) %>%
    dplyr::filter(df_counts_mtx$Sample.ID == sample_id) %>%
    group_by(fov) %>%
    summarize(count = n(), # number of cells
              Cdkn1a_count = sum(Cdkn1a > 0),
              Mki67_count = sum(Mki67 > 0),
              Lmna_count = sum(Lmna > 0), 
              .groups = "drop") %>%
    dplyr::filter(count > 0)  %>%
    mutate(sample=sample_id,
           pct_cells_Cdkn1a = (Cdkn1a_count / sum(count)) * 100,
           pct_cells_Mki67 = (Mki67_count / sum(count)) * 100,
           pct_cells_Lmna = (Lmna_count / sum(count)) * 100) %>%
    dplyr::relocate(sample) 
  
  return(df_pct_fov)
}



pct_non_zero <- function(x) {
  length(x[x >0]) / length(x) * 100
}

mean_non_zero <- function(x) {
  mean(x[x >0])
}

get_top_markers <- function(sobj, cluster, n=10) {
  markers <- FindMarkers(sobj, ident.1 = cluster)
  markers %>%
    rownames_to_column(var = "gene") %>% # Add gene names as a column
    dplyr::top_n(n, wt = abs(avg_log2FC)) %>%  # Select top 10 markers by absolute log2 fold change
    dplyr::mutate(cluster = cluster)            # Add cluster information
}



# DotPlot
my_dot_plot <- function (sobj, cluster_col='seurat_clusters', title='', results_dir, genes_cosMx_celltyping) {
  
  # Function to find top n markers for a specific cluster
  clusters <- unique(sobj@meta.data[[cluster_col]])
  print(clusters)
  # Apply the function to all clusters and combine results
  Idents(sobj) <- cluster_col
  top_markers <- map_dfr(clusters, get_top_markers, sobj=sobj, n=15)
  
  #table(top_markers$cluster)
  dim(top_markers) # num clusters * n
  print(head(top_markers))
  # also append our favorite cell typing markers
  
  # Genes we want to see because they are classic cell typing genes
  #genes_lineage <- c("Cd19", "Ptprc", "Cd8a", "Cd8b1", "Cd68", "Cd163","Itgax", "Itgam", "Cd3e", "Cd3d", "Fn1", "Acta2")
  
  # Compare DotPlots with top_markers limited by only ns celltyping genes, but
  # augmented with classic cell typing genes for comparison
  # remove any genes where the rowsum is too low unless its a lineage marker, try to keep with lower barriers. 
  markers_to_plot <- top_markers %>%
    dplyr::filter ((top_markers$gene %in% genes_cosMx_celltyping & (top_markers$pct.1 + top_markers$pct.2) > 0.1
                    & top_markers$p_val_adj < 0.05) | top_markers$gene %in% genes_lineage)

  dim(markers_to_plot) # 72 is somewhat reasonable; w repeats is really 36 unique genes

  p <- DotPlot(sobj, features = unique(markers_to_plot$gene)) + RotatedAxis() +
    theme_minimal() + # Change to a minimal (white) background
    theme(
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1), # Reduce font size and rotate x-axis labels
      axis.text.y = element_text(size = 10), # Adjust y-axis font size if needed
      panel.background = element_rect(fill = "white", color = NA), # Force white panel background
      plot.background = element_rect(fill = "white", color = NA)  # Force white overall background
      #panel.grid = element_blank() # Optional: remove grid lines if needed
    ) +
    ggtitle(title)

  ggsave(paste0(results_dir, "/",sanitize_name(title),"_DotPlot.png"), p, width = 12)
  return(p)
}

# Get Percentage of cells for each fov that are positive, independent of celltype
get_counts_byfov <- function(df, count_matrix, df_rna_lookup, codeclass){
  
  # Calculate the total counts per cell (column sums)
  filtered_display_names <- df_rna_lookup$DisplayName[df_rna_lookup$CodeClass == codeclass]
  cat(paste("number of genes to include:",length(filtered_display_names)))
        
  filtered_counts_matrix <- counts_matrix[rownames(counts_matrix) %in% filtered_display_names, ]
  dim(filtered_counts_matrix)  
  
  total_counts <- colSums(filtered_counts_matrix)
  
  # Convert to a dataframe
  total_counts_df <- data.frame(
    cell_id = names(total_counts),
    num_probes = total_counts
  )
  # Set the rownames to the cell IDs
  rownames(total_counts_df) <- total_counts_df$cell_id
  
  # merge with df
  df2 <- merge(df, total_counts_df, by='cell_id')
  
  df_counts <- df2 %>%
    dplyr::select(c('Sample.ID', 'fov', 'num_probes')) %>%
    group_by(fov, Sample.ID) %>%
    summarize(cell_count = n(), # number of cells
              num_probes = sum(num_probes)) # %>%
              #.groups = "drop" 
   # dplyr::relocate(sample) 
  
  return(df_counts)
}



