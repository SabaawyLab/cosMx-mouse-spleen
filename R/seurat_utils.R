
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

# TODO: Needed?

#'  write out the expression matrix from a Seurat object to file
#'
#' @param seurat_obj The Seurat object
#' @param filename The filename for the matrix data
#' @return success boolean
#' @examples
#' success <- export_seurat_expr_matrix(seurat_obj, "data/seurat_counts.csv")
# export_seurat_expr_matrix <- function(seurat_obj, filename) {
#   # write out the expression matrix
#   # TODO: use the Default Assay rather than hardcoding RNA
#   write.csv(as.matrix(seurat_obj@assays$RNA@counts), file = filename)
#   return(true)
# }

# export_seurat_metadata <- function(seurat_obj, filename) {
#   # write out the metadata
#   metadata <- get_metadata(seurat_obj)
#   write.csv(metadata, file = filename)
# }
# 
# export_seurat_spatial <- function(seurat_obj, filename) {
#   # write out the spatial data
#   spatial <- seurat_obj@spatial
#   write.csv(spatial, file = filename)
# }


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


flip_centroids <- function(df_meta){
  min_y <- min(df_meta$CenterY_global_px)
  max_y <- max(df_meta$CenterY_global_px)
  df_meta$CenterY_global_px <- (max_y - df_meta$CenterY_global_px) + min_y
  return(df_meta)
}

#' Load metadata flat file from cosMx flatfiles, 
#' given Seurat LoadNanostring() doesn't load this
#' Assumptions about columns existing in the nanostring format, e.g. fov, cell_ID exist
#' Rows do not need to match precisely, but cell_ID must match as an integer.
#'
#' @param sobj The Seurat object
#' @param metadata_file The filename for the metadata file
#' @return The Seurat object with metadata loaded
#' @examples
#' sobj <- load_meta_load_missed(seurat_obj, "flatfiles/metadata.csv")
load_meta_load_missed <- function(sobj, metadata_file) {
  # Seurat LoadNanostring() workaround to load metadata flat file including making sure the keys match
  # this loads morph markers, measurements, X, Y, (global and local) and cell_ID
  # TODO: add error checking on cell_ID format as integer
  
  #df_obj6_meta <- read.csv(paste0(lung6.data.dir, "/Lung6_metadata_file.csv"))
  df_obj_meta <- read.csv(metadata_file)
  
  # workaround to fix inverted centroids, flips CenterY_global_px
  #df_obj_meta <- flip_centroids(df_obj_meta)
  
  # shorten unweildy column names
  string_to_replace <- "fef9b11c.70b7.4905.b96d.0b80c5560fb7"
  replacement_string <- "fef"
  names(df_obj_meta) <- gsub(string_to_replace, replacement_string, names(df_obj_meta))
  string_to_replace <- "61b859e8.677d.4cad.9cf5.44e03fe9960e"
  replacement_string <- "61b"
  names(df_obj_meta) <- gsub(string_to_replace, replacement_string, names(df_obj_meta))
  
  df_obj_meta$cell_num = df_obj_meta$cell_ID # save integer version of cell_ID
  #df_obj_meta$key = paste0("c_", slide_num, "_",df_obj_meta$fov, "_", df_obj_meta$cell_ID)
  df_obj_meta$key = paste0( df_obj_meta$cell_ID, "_", df_obj_meta$fov)
  rownames(df_obj_meta) <- df_obj_meta$key
  
  head(rownames(sobj@meta.data))
  print("Number of cells in Seurat object before adding metadata")
  print(nrow(sobj@meta.data))
  # AddMetaData needs the same records in the same order
  print("Do the number of rows in the metadata file match the sobj?")
  print(nrow(df_obj_meta) == nrow(sobj@meta.data)) # FALSE
  lost_on_load <- setdiff(rownames(df_obj_meta), rownames(sobj@meta.data)) # about 100, ok, lets drop these 
  length(lost_on_load)
  if(length(lost_on_load) > 0){
    print ("These are the rows that are in the metadata file but not the Seurat obj, we will drop them from the metadata file.")
    print(head(lost_on_load))
  }
  print("expect 0 cell differences")
  print(setdiff(rownames(sobj@meta.data), rownames(df_obj_meta))) # 0
  
  # delete these from df_obj_metadata
  df_obj_meta <- df_obj_meta[!rownames(df_obj_meta) %in% lost_on_load,]
  print(setdiff(rownames(df_obj_meta), rownames(sobj@meta.data))) # 0 now.
  
  df_obj_meta$cell_ID <- df_obj_meta$key # restore cell_ID in standard format
  sobj <- AddMetaData(sobj, metadata = df_obj_meta)
  head(sobj@meta.data)
  head(rownames(sobj@meta.data))
  
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
      print("The keys match but are not in the same order")
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


#' check the Seurat object for size, health, and status;
#' report on the number of cells, genes, and metadata
#'
#' @param sobj The Seurat object
#' @return success boolean, False if something is very wrong; 0 cells or genes, or NAs in rownames
#' @examples
#' success <- check_seurat_obj(sobj)
check_seurat_obj <- function(sobj){
  # check the Seurat object for size, health, and status
  # report on the number of cells, genes, and metadata
  # TODO: add more checks relating to cosMx data expectations
  print(paste("Version", sobj@version))
  print(paste("Number of cells", nrow(sobj@meta.data)))
  print("Expect no NAs in the rownames:")
  print(sum(is.na(rownames(sobj@meta.data))))
  print(paste("Layers", Layers(sobj)))
  print(paste("Assays", Assays(sobj)))
        
  print(paste("Number of genes", length(Features(sobj))))
  
  print("Number of metadata columns")
  print(ncol(sobj@meta.data))
  
  # What FOVs are defined
  print(sobj@images)
  
  # Layer for counts
  # Assay for Nanostring

  # Do the important metadata columns exist? rownames, cell_id, global X & Y coords, local X & Y coords, fov
  
  # Do we have an FOV defined that matches the metadata? This will enable ImageDimPlot() to work. 
  
  # anything NA that shouldn't be?  
  
  # expect 960 features
  print(paste0("Num Features expect 960: ", nrow(sobj)))
  
  return (TRUE)
}


#' subtract the negative controls from the counts matrix
#' and return the modified sparse count matrix, 0 values in the Negative features 
#' Assumes the Negative features are still present in the Seurat object 
#' and the Assay is "Nanostring", layer is "counts", per defaults from Seurat::LoadNanostring()
#'
#' @param obj The Seurat object
#' @return modified count matrix
#' @examples
#' obj <- subtract_neg_control(obj)
subtract_neg_control <- function(obj){
  # get negative probe feature indexes, pre-normalization 
  counts_matrix <- GetAssayData(object = obj, assay = "Nanostring", layer="counts") ## or should use normalized?
  
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
  obj <- SetAssayData(object = obj, assay = "Nanostring", layer="counts", new.data = count_replace)
  
  return(obj)

}

#' remove the system control features from the counts matrix
#' and return the modified Seurat object.  
#' Tip: Do this before SCTransform to avoid needing to subtract these from more layers
#'
#' @param obj The Seurat object
#' @param assay The assay to use, default is "Nanostring"
#' @param layer The layer to use, default is "counts"
#' @return sobj with system control features removed
#' @examples
#' obj <- remove_sys_control(obj)
#' obj <- remove_sys_control(obj, assay="Nanostring", layer="counts")
remove_sys_control <- function(obj){
  # get negative probe feature indexes, pre-normalization 
  counts_matrix <- GetAssayData(object = obj, assay = "Nanostring", layer="counts") ## or should use normalized?
  
  feature_names <- rownames(obj)
  # which feature_names in list begin with "Neg"?
  
  feature_names 
  sys_control_names <- feature_names[grep("^SystemControl", feature_names)]
  genes_to_keep <- feature_names[!feature_names %in% sys_control_names]
  
  obj <- subset(obj, features = genes_to_keep)

  # Get row indices 
  #keep_indices <- which(!(feature_names %in% sys_control_names))
  # Create sparse matrices for neg_feature_names and others
  #keep_matrix <- counts_matrix[keep_indices, ]
  
  #print("Number of system control features to remove:")
  #print(length(sys_control_names))
  print(paste("Remaining dimensions:", paste(dim(obj), collapse = " x ")))
  #print(dim(keep_matrix))   # Dimensions of the matrix without 
 
  #obj[[assay]]@counts <- keep_matrix
  
  return(obj)
  
}

#' retrieve the cell-wise mean of the negative control features
#' Assumes: Negative control features are still present in the Seurat object
#'
#' @param obj The Seurat object
#' @return numeric vector of mean negative control values by cell_id
#' @examples
#' obj <- get_neg_control_means(obj)
get_neg_control_means <- function(obj){
  
  counts_matrix <- GetAssayData(object = obj, assay = "Nanostring", layer="counts") ## or should use normalized?
  
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

## Temporary PanCK/CD45 handling
# Define the custom function
# Note: These hard coded values were obtained by manually observing the Mean intensity values in the component images
calculate_PanCK_PT <- function(value) {
  # if (is.na(value) || value > 4000) {
  #   return('na')
  if (value > 700) {
    return(1)
  } else {
    return(0)
  }
}

calculate_CD45_PT <- function(value) {
  if (is.na(value)) {
    return(0)
  } else
    if (value > 200) {
      return(1)
    } else {
      return(0)
    }
}


#' Add FOV and Sample metadata to the Seurat object
#'
#' @param obj The Seurat object
#' @df_fov_meta dataframe with FOV and Sample metadata with columns 'fov' and 'Run_Tissue_name' to match the object's metadata as keys
#' @return Seurat Object
#' @examples
#' obj <- add_sample_metadata(obj, df_fov_meta)
add_sample_metadata <- function(obj, df_fov_meta){
  
  df_metadata <- obj@meta.data
  # Apply the function to create the PanCK.PT column
  
  
  df_metadata2 <- merge(df_metadata, df_fov_meta, 
                        by.x = c('fov', 'Run_Tissue_name'), by.y = c('fov', 'Run_Tissue_name'),
                        suffixes = c(".x",""), how = 'inner')
  #print(head(df_metadata2))
  
  df_metadata2$Mean.PanCK[is.na(df_metadata2$Mean.PanCK)] <- 0
  df_metadata2$Max.PanCK[is.na(df_metadata2$Max.PanCK)] <- 0
  #df_metadata <- df_metadata %>%
  #  mutate(PanCK.PT = sapply(Mean.PanCK, calculate_PanCK_PT))
  
  # sample thresholds
  panCK.thresholds <- c("1_1" = 5000, "2_4" = 3000, "3_5" = 2000, "3_6" = 4000, "4_1" = 1800, "4_6" = 1300,
                        "5_2" = 1000, "6_3" = 2000, "6_5" = 10000, "7_4" = 10000,
                        "ORG_2_3" = 0, "ORG_4_1" = 0, "ORG_8_2" = 0, "ORG_8_5" = 0,
                        "A1" = 0, "A2" = 0, "A3" = 0, "B1" = 0, "B2" = 0, "B3" = 0, 
                        "C1" = 0, "C2" = 0)
  df_metadata2 <- df_metadata2 %>%
    dplyr::mutate(
      PanCK.threshold = panCK.thresholds[Sample.ID],  # Look up threshold by Sample.ID
      PanCK.PT = as.integer(Mean.PanCK >= PanCK.threshold)  # Threshold comparison
    )
  
  df_metadata2$Mean.CD45[is.na(df_metadata2$Mean.CD45)] <- 0
  df_metadata2$Max.CD45[is.na(df_metadata2$Max.CD45)] <- 0
  CD45.thresholds <- c("1_1" = 300, "2_4" = 200, "3_5" = 200, "3_6" = 200, "4_1" = 200, "4_6" = 200,
                       "5_2" = 200, "6_3" = 1000, "6_5" = 200, "7_4" = 230,
                       "ORG_2_3" = 1, "ORG_4_1" = 1, "ORG_8_2" = 1, "ORG_8_5" = 1,
                       "A1" = 1, "A2" = 1, "A3" = 1, "B1" = 1, "B2" = 1, "B3" = 1, "C1" = 1, "C2" = 1)
  
  #df_metadata <- df_metadata %>%
  #  mutate(CD45.PT = sapply(Mean.CD45, calculate_CD45_PT))
  df_metadata2 <- df_metadata2 %>%
    mutate(CD45.threshold = CD45.thresholds[Sample.ID],          # Look up threshold by sample
           CD45.PT = as.integer(Mean.CD45 >= CD45.threshold)) 
  
  cols_to_add <- c("cell_id", "CD45.PT", "PanCK.PT", "Mean.PanCK", "Mean.CD45", "Diagnosis",
                   "PanCK.threshold", "CD45.threshold", colnames(df_fov_meta))
  print(paste("cells before merge", nrow(df_metadata2)))
  
  df_metadata2$PanCK.PT[is.na(df_metadata2$PanCK.PT)] <- 1
  df_metadata2$CD45.PT[is.na(df_metadata2$CD45.PT)] <- 0

  print(paste("cells after merge", nrow(df_metadata2)))
  
  rownames(df_metadata2) <- df_metadata2$cell_id 
  # keep only some columns
  df_metadata2 <- df_metadata2[, cols_to_add]
  
  # Do we have all the slides?
  table(df_metadata2$Run_Tissue_name)
  nrow(df_metadata2) # same as full orig. ?
  
  # Identify rows with NA in cell_ID
  na_rows <- is.na(df_metadata2$cell_ID)
  # Display rows with NA in cell_ID
  nrow(df_metadata2[na_rows, ])
  
  # do_rownames_match(sobj, df_metadata2) # FALSE 
  # display cell_Ids in sobj not found in df_metadata2
  cell_IDs_not_in_metadata <- setdiff(colnames(obj), rownames(df_metadata2))
  length(cell_IDs_not_in_metadata) # expect 0 
  #cell_IDs_not_in_metadata[1:10]
  cell_IDs_not_in_sobj <- setdiff(rownames(df_metadata2), colnames(obj)) # 0
  # sort 
  df_metadata2 <- df_metadata2[colnames(obj), ]
  do_rownames_match(obj, df_metadata2) # try again?
  
  # preview metadata for cell_ID "c_4_18_1"
  #df_metadata[df_metadata$cell_ID == "c_4_18_1",]
  
  obj <- AddMetaData(obj, df_metadata2)
  
  return(obj)
}


#' Combine two Cell Profile Matrices for use with InSituType
#' Assumes the genes are the same in both matrices, and the Cell Types are Mutually Exclusive
#' For example, use this to append specific tumor subtypes to a general healthy cell type matrix
#'
#' @param m1 a count matrix with cells as colnames and genes as rownames
#' @param m2 a matrix with the same genes as m1
#' @return numeric vector of mean negative control values by cell_id
#' @examples
#' new_matrix <- merge_matrices(matrix1, matrix2)
merge_matrices <- function(m1, m2){
  # Given two InsituType matrices, merge them by rownames
  # TODO: add test method
  df1 <- as.data.frame(m1)
  df2 <- as.data.frame(m2)
  df1$genes <- rownames(df1)
  df2$genes <- rownames(df2)
  merged <- merge(df1, df2, by = "genes", all = TRUE)
  rownames(merged) <- merged$genes
  merged <- merged[, -1]
  merged[is.na(merged)] <- 0
  
  # insitu_matrix <- as.matrix(insitu_matrix_He)
  # rownames(insitu_matrix) <- merged$genes
  # insitu_matrix <- apply(insitu_matrix, 2, as.numeric)
  # rownames(insitu_matrix) <- rownames(merged)
  # rownames(insitu_matrix)[1:10]
  # head(insitu_matrix)
  # pretty sure the NAs are for genes we won't have anyway, but..
  return(merged)
}


#' Grouped plot of cell counts, colored by cell type 
#'
#' @param cell_counts an aggregated data frame of cell counts by 2 groupings
#' @param group.by the column name to group by
#' @param title a string to append to the Title of the plot
#' @return numeric vector of mean negative control values by cell_id
#' @examples
#'  plot_celltype(grouped_df)
plot_celltypes <- function(cell_counts_df, group.by="Sample.ID", title="by Sample.ID",colors = c(""), log2=FALSE){
  
  # cell_counts <- df %>%
  #   group_by(Run_Tissue_name, cell_type) %>%
  #   summarise(count = n(), .groups = 'drop')
  print(head(cell_counts_df))
  
  num_levels <- length(unique(cell_counts_df$cell_type))
  #dynamic_colors <- brewer.pal(min(8, num_levels), "Dark2")  # Adjust as needed, max typically 8-12 for most palettes
  # Or use rainbow? (hideous)
  #dynamic_colors <- rainbow(num_levels)
  
  #colors() has many nice ones. 
  if (length(colors) == 1) {
    colors <- c("darkblue", "red", "orange", "#996633" ,"green", 
                 '#B3DE69', "#FFFF00", 
                 "#FF6633", '#FDB462',
                 "yellowgreen", "springgreen", "#006600" ,
                 "#999999" ,
                 "#33CC00"  ,"#9966CC"  ,"grey10", "#3399CC"  ,"#FF0000"  , "#CC0000", "lightblue",
                 "#FF9900" ,  "#FF66FF" , 'purple', '#8DD3C7', '#FB8072', '#80B1D3',
                  'magenta', 'blue', '#BC80BD','cyan', 'green',
                 "turquoise1", "red", "yellow",  "green", "black", "grey50", "thistle", 
                 "salmon3", "tan", "darkorange","wheat", "violet",  "rosybrown2")
  
    names(colors) <- c("endothelial", "tumor", "epithelial", "fibroblast", "immune",
                        "lymphocyte", "macrophage", 
                        "adeno", "SqCC",
                        "innate immune", "T.cell",  
                        "B.cell", "other", 
                        "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", 
                        "t", "u", "v", "w", "x", "y", "z","Other", "Annotated Tumor")
  }
  if (log2) {
    cell_counts_df$count <- log2(cell_counts_df$count + 1)
  }
  # Create a bar plot
  p <- ggplot(cell_counts_df, aes(x = .data[[group.by]], y = count, fill = cell_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Run Tissue Name", y = "Count of Cells") + # , fill = "Cell Type") +
    #scale_fill_brewer(palette = "Dark2") + # also try Set2
    # another option: 
    scale_fill_manual(values = colors)  # Adjust colors as needed based on the number of cell types
    #theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))   # Rotate x-axis labels for better readability
  p <- p + ggtitle(paste("Cell Type Counts by", group.by, title))
  
  # Display the plot
  print(p)
  
}

#' Given a seurat object (e.g. maybe a subset), a vector of negative control means, and a title suffix,
#' Run InSituType and plot the results
#' Returns the cell_type assignments by cell_id as a "sup" list, aka, what InSituType returns
#'
#' @param obj 
#' @param negmeans a labelled list of negative control means by cell_id
#' @param title_suffix a string to append to the Title of the plot
#' @param group.by a string to group the cells by e.g. Run_Tissue_name
#' @param insitu_matrix a matrix of Cell Type average gene expressions: "Cell Profile Matrix"
#' @param results_dir the path to save plots
#' @return sup a list of cell types by cell id
#' @examples
#' run_intitu_type (sobj_fov24, negmenas, "lung FOV24 only", "cell_type", insitu_matrix_He)
run_insitu_type <- function(obj, negmeans, title_suffix="", group.by="Sample.ID", insitu_matrix, results_dir){
  
  # align cells in negmeans with cells in obj
  df_meta = obj@meta.data
  # Get the mean of the negative controls, already calculated above. 
  
  counts_matrix <- GetAssayData(object = obj, assay = "SCT", layer="counts")
  counts <- counts_matrix %>%
    as.matrix() %>%
    t()
  

  if (sum(negmeans) == 0) {
    stop("Negative control means are all zero. Please check the negative control names.")
  }
  negmeans[1:10]
  # only keep cells that match the rownames of counts matrix 
  negmeans <- negmeans[rownames(counts)]
  
  if (group.by == ""){
    group.by <- "Run_Tissue_name"
  }
  # Run InSituType
  sup <- insitutypeML(x = counts,
                      neg = negmeans,
                      cohort = NULL, # cohort_data,
                      reference_profiles = insitu_matrix)  
  
  # Plot the results
  # Plot the AtoMx Cell type breakdowns by slide
  # "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"  # first one, supervised, about 30 types
  # "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" # unsupervised, 20 clust
  # RNA_nbclust_c5e1d3d8.0a17.443a.994b.9a20f297fd2c_1_clusters"   
  df_meta$cell_type <- sup$clus
  
  df_meta <- df_meta %>%
    mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
    select(!!group.by, cell_type)
  
  # df_meta <- df_meta %>%
  #   mutate(Run_Tissue_name = gsub(" ", "_", Run_Tissue_name)) %>%
  #   select(Run_Tissue_name, cell_type)
  # 
  df_celltype <- df_meta %>%
    group_by_at(c(group.by, "cell_type")) %>%
    summarise(count = n(), .groups = 'drop')
  
  p <- plot_celltypes(df_celltype, group.by, title_suffix)
  ggsave(paste(results_dir, paste0("CellType_by_",group.by,"_",sanitize_name(title_suffix),".png")), plot=p, width = 10, height = 6, dpi = 300)
  
  # requires PCA to be run
  #p2 <- DoHeatmap(obj, features=TopFeatures(obj, dim=1, nfeatures=20))
  #ggsave(paste(results_dir, paste0("InSitu_Heatmap_for_",sanitize_name(title_suffix),".png")), 
  #       plot=p2, width = 10, height = 6, dpi = 300)
  
  return(sup)
  
}
  

run_insitu_type2 <- function(obj, negmeans, n_clusts = 0,title_suffix="", group.by="Sample.ID", insitu_matrix, results_dir){

  set.seed(0)
  # align cells in negmeans with cells in obj
  df_meta = obj@meta.data
  # Get the mean of the negative controls, already calculated above. 
  
  counts_matrix <- GetAssayData(object = obj, assay = "SCT", layer="counts")
  counts <- counts_matrix %>%
    as.matrix() %>%
    t()
  
  
  if (sum(negmeans) == 0) {
    stop("Negative control means are all zero. Please check the negative control names.")
  }
  negmeans[1:10]
  # only keep cells that match the rownames of counts matrix 
  negmeans <- negmeans[rownames(counts)]
  
  if (group.by == ""){
    group.by <- "Run_Tissue_name"
  }
  # Run InSituType
  sup <- insitutype(x = counts,
                      neg = negmeans,
                      n_clusts = n_clusts,
                      cohort = NULL, # cohort_data,
                      reference_profiles = insitu_matrix,
                      update_reference_profiles = TRUE
                    )  
  
  # Plot the results
  # Plot the AtoMx Cell type breakdowns by slide
  # "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"  # first one, supervised, about 30 types
  # "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" # unsupervised, 20 clust
  # RNA_nbclust_c5e1d3d8.0a17.443a.994b.9a20f297fd2c_1_clusters"   
  df_meta$cell_type <- sup$clus
  
  df_meta <- df_meta %>%
    mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
    select(!!group.by, cell_type)
  
  # df_meta <- df_meta %>%
  #   mutate(Run_Tissue_name = gsub(" ", "_", Run_Tissue_name)) %>%
  #   select(Run_Tissue_name, cell_type)
  # 
  df_celltype <- df_meta %>%
    group_by_at(c(group.by, "cell_type")) %>%
    summarise(count = n(), .groups = 'drop')
  
  p <- plot_celltypes(df_celltype, group.by, title_suffix)
  ggsave(paste(results_dir, paste0("CellType_by_",group.by,"_",sanitize_name(title_suffix),".png")), plot=p, width = 10, height = 6, dpi = 300)
  
  # requires PCA to be run
  #p2 <- DoHeatmap(obj, features=TopFeatures(obj, dim=1, nfeatures=20))
  #ggsave(paste(results_dir, paste0("InSitu_Heatmap_for_",sanitize_name(title_suffix),".png")), 
  #       plot=p2, width = 10, height = 6, dpi = 300)
  
  return(sup)
  
}

#' Given a seurat object (e.g. maybe a subset), the semi-supervised or supervised cell type assignments,
#' and the results directory
#' Run InSituType and plot the results
#' Returns the cell_type assignments by cell_id as a "sup" list, aka, what InSituType returns
#'
#' @param obj 
#' @param insitutype_results a labelled list of negative control means by cell_id, needed for flightpath
#' @param results_dir the path to save plots
#' @param colorset a named (or unnamed) vector of colors for each cell type
#' @param group.by a string to group the cells by e.g. Run_Tissue_name
#' @param title_suffix a string to append to the Title of the plot
#' @param cell_labels_col name of the column in metadata with cell type labels

#' @return sup a list of cell types by cell id
#' @examples
#' run_intitu_type (sobj_fov24, negmenas, "lung FOV24 only", "cell_type", insitu_matrix_He)
#' # TODO: use results object for plotting fp
insitu_cell_typing_plots <- function(sobj, 
                                     insitutype_results, 
                                     results_dir, 
                                     colorset,
                                     group.by = "Run_Tissue_name", 
                                     title_suffix = "InSituType",
                                     cell_label_col = "CellType_Labels" 
                                     ){
  
  # does metadata column exist?
  if (!cell_label_col %in% colnames(sobj@meta.data)){
    stop("The column ", cell_label_col, " does not exist in the metadata")
  }
  
  if (missing(colorset)){
     cols <- InSituType::colorCellTypes(freqs = table(sobj@meta.data[[cell_label_col]]), palette = "brewers")
  } else {
    cols <- colorset
  }
  # if no attempt at naming was made, fill in
  if (is.null(names(colorset))){
    print("Assinging colors sequentially")
    names(cols) <- unique(sobj@meta.data[[cell_label_col]])
    #names(cols) <- unique(insitutype_results$clust)
  }
  # if some items only remain unnamed, fill in with defaults from brewers
  cols1 <- InSituType::colorCellTypes(freqs = table(sobj@meta.data[[cell_label_col]]), palette = "brewers")
  cols <- cols[seq_along(unique(sobj@meta.data[[cell_label_col]]))]
  names(cols) <- unique(sobj@meta.data[[cell_label_col]])
  # 2. match up named colors in iocolors + colors_tumor1
  if (!missing(colorset)){
    cols[is.element(names(cols), names(colorset))] <- colorset[names(cols)[is.element(names(cols), names(colorset))]]
  }

  # are all labels in the colorset?
  # labels_missed <- setdiff(unique(sobj@meta.data[[cell_label_col]]), names(cols))
  # if (length(labels_missed) > 0){
  #   print(c("labels needed in metadata: ", labels_missed))
  #   stop("Not all cell types are in the colorset")
  # }
  print(cols)
  
  # #cols <- colorset
  # str(cols)
  # colors <- colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(21)
  
  # This plot needs a hardcoded column named cell_type
  sobj <- AddMetaData(sobj, 
                       metadata = data.frame(cell_type = sobj@meta.data[[cell_label_col]],
                                             cell_id = rownames(sobj@meta.data)))
  df_meta <- sobj@meta.data
  
  print(table(df_meta$cell_type))
  
  df_meta <- df_meta %>%
    mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
    select(!!group.by, cell_type)
  
  df_celltype <- df_meta %>%
    group_by_at(c(group.by, "cell_type")) %>%
    summarise(count = n(), .groups = 'drop')
  
  p <- plot_celltypes(df_celltype, group.by, col=colorset, title_suffix)
  ggsave(paste0(results_dir, paste0("/CellType_by_",group.by,"_",sanitize_name(title_suffix),".png")), plot=p, width = 10, height = 6, dpi = 300)
  
  # flightpath
  

  #>  Named chr [1:21] "#8DD3C7" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" ...
  #>  - attr(*, "names")= chr [1:21] "f" "b" "a" "c" ...
  print(cols[c("Endothelial", "Fibroblast", "B.cell")])
  
  # flightpath <- InSituType::flightpath_layout(logliks = insitutype_results$logliks, 
  #                                             profiles = insitutype_results$profiles)
  # 
  # # Save the plot as a PNG file
  # png(paste0(results_dir, paste0("/InSitu_FlightPath_for_",sanitize_name(title_suffix),".png")), width = 600, height = 600)  # Open the PNG device
  # par(mar = c(0,0,0,0))
  # plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[insitutype_results$clust])
  # text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
  # dev.off()
  #ggsave(paste0(results_dir, paste0("/InSitu_FlightPath_for_",sanitize_name(title_suffix),".png")), width = 10, height = 6, dpi = 300)
  
  #this assumes sobj has umap coordinates
  # plot(um, pch = 16, cex = 0.1, col = cols[sup$clust])
  # for (cell in unique(insitutype_results$clust)) {
  #   text(median(um[insitutype_results$clust == cell, 1]), median(um[insitutype_results$clust == cell, 2]), cell, cex = 0.7)
  # }
  
  # spatial
  # xy space, by slide or sample
  # how to group semisup10 by sample
  sample.groups <- setNames(sobj@meta.data$Sample.ID, sobj@meta.data$cell_id)
  sample.groups[1:5]
  head(insitutype_results$clust)
  combined_df <- data.frame(insitutype_results$clust, sample.groups, 
                            row.names = names(insitutype_results$clust))
  head(combined_df)
  

  # for each sample.group in combined_df
  #   plot the xy points for that sample.group
  #   color the points by semisup10$clust
  
  # Assuming sample.groups is a vector and combined_df is your dataframe
  unique_groups <- unique(sample.groups)
  
  df_meta <- sobj@meta.data
  # HARDcoded
  combined_df <- df_meta[c('CellType_Labeled', 'Sample.ID')]
  colnames(combined_df) <- c('semisup.clust', 'sample.groups')
  # Loop over each unique group
  for (group in unique_groups) {
    # Subset the combined_df for rows where sample.groups matches the current group
    group_subset <- combined_df[sample.groups == group, ]
    
    # Perform operations with group_subset, e.g., print the subset
    print(paste("Group:", group))
    print(paste("group_subset cells ", nrow(group_subset)))
    #plot(xy, pch = 16, cex = 0.1, col = cols[group_subset$semisup10.clust], asp = 1) 
    
    group_subset <- combined_df[sample.groups == group, ]
    # Perform operations with group_subset, e.g., print the subset
    print(paste("Group:", group))
    print(paste("group_subset cells ", nrow(group_subset)))
    df_meta_grp <- df_meta[df_meta$Sample.ID == group, ]
    
    # filter down df_meta_grp to exclude cell types with < 10 cells
    df_celltype <- df_meta_grp %>%
      group_by(CellType_Labeled) %>%
      summarise(count = n(), .groups = 'drop') %>%
      mutate(percentage = (count / sum(count)) * 100)
    
    # filter down df_meta_grp to exclude cell types with < 10 cells
    df_meta_grp <- df_meta_grp %>%
      filter(CellType_Labeled %in% df_celltype$CellType_Labeled[df_celltype$count >= 20])
    
    p <- ggplot(df_meta_grp,
                aes(CenterX_global_px, y=CenterY_global_px, color=cell_type)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(values = cols[df_meta$cell_type]) +
      coord_fixed(ratio = 1) +
      ggtitle(paste("Sample", group))
    ggsave(paste0(results_dir, "/sample_", group,"_",sanitize_name(title_suffix),"_xy.png"), p)
  }

  return(TRUE)
  
}


plot_morphology_markers <- function(df_meta, 
                                    group.by = "Run_Tissue_name", 
                                    title_suffix = "Morphology Mkrs") {
  
  unique_groups <- unique(df_metadata[[group.by]])
  
  for (group in unique_groups) {
    # Subset the combined_df for rows where sample.groups matches the current group
    group_subset <- df_meta[df_meta[[group.by]] == group, ]
    
    # Perform operations with group_subset, e.g., print the subset
    print(paste("Group:", group))
    print(paste("group_subset cells ", nrow(group_subset)))
    #plot(xy, pch = 16, cex = 0.1, col = cols[group_subset$semisup10.clust], asp = 1) 
    
    group_subset <- df_meta[df_meta[[group.by]] == group, ]
    # Perform operations with group_subset, e.g., print the subset
    print(paste("Group:", group))
    print(paste("group_subset cells ", nrow(group_subset)))
    df_meta_grp <- df_meta[df_meta[[group.by]] == group, ]
    
    # set color to green and magenta, w legen
    p <- ggplot(df_meta_grp,
                aes(CenterX_global_px, y=CenterY_global_px, color=cell_type)) +
      geom_point() +
      scale_color_manual(values = cols[df_meta$cell_type]) +
      coord_fixed(ratio = 1) +
      ggtitle(paste("Sample", group))
  }
}

# Function to calculate the average by group
#'
#' @param matrix string with gene name
#' @param cell_type_list 
#' @return named vector of separated gene names
#' @examples
#' group_averages <- calculate_group_averages(non_neg_matrix, cell_type_list)
calculate_group_averages <- function(matrix, groups) {
  unique_groups <- unique(groups)
  group_averages <- sapply(unique_groups, function(group) {
    group_indices <- which(groups == group)
    group_matrix <- matrix[, group_indices, drop = FALSE]
    rowMeans(group_matrix)
  })
  return(group_averages)
}




#' Reformat gene names that have been combined with a slash, 
#' e.g. CXCL1/2/3 is 3 genes: CXCL1, CXCL2, CXCL3
#' EIF5A/L1 is 2 genes: EIF5A, EIF5AL1 but we'll find only EIF5A (ok)
#'
#' @param gene string with gene name
#' @return named vector of separated gene names
#' @examples
#' expanded_genes <- expand_gene_names("CXCL1/2/3")
expand_gene_names <- function(gene) {
  # first, remove dash or spaces
  gene <- gsub("-", "", gene)
  gene <- gsub(" ", "", gene)
  
  # Split the gene name by '/', if any
  parts <- str_split(gene, "/", simplify = TRUE)
  
  if (length(parts) == 1) {
    return(parts)
  }
  #print(parts)
  # check pattern one or two characters per part
  #print(nchar(parts[2]))
  if (nchar(parts[2]) > 1) {
    #change parts[1] to have same length as parts[2], taken from end of parts[1] string
    #parts[1] <- substr(parts[1], 1, nchar(parts[2]))
    first_div <- substr(parts[1], nchar(parts[1]) - 1, nchar(parts[1]))
    prefix <- substr(parts[1], 1, nchar(parts[1]) - 2)
    #print(paste("first_div", first_div, "prefix", prefix))
    parts[1] <- first_div
  }  
  else if (nchar(parts[2]) == 1) {
    first_div <- substr(parts[1], nchar(parts[1]), nchar(parts[1]))
    prefix <- substr(parts[1], 1, nchar(parts[1]) - 1)
    #print(paste("first_div", first_div, "prefix", prefix))
    parts[1] <- first_div
  }
  # expect to now have gene_part and parts as vector of suffixes e.g. A, B, C or L1, L2
  
  # Get the prefix before the last occurrence of non-slash characters before a slash
  #print(prefix)
  # Reconstruct each part with the prefix
  expanded_names <- paste0(prefix, parts)
  # Return a single string with each part separated by commas
  return(paste(expanded_names, collapse = ", "))
}

# helper function to reformat M:M vector of gene equivalencies to a dataframe
gene_split_to_df <- function(expanded_gene_vector){
  
  # Convert the named vector to a dataframe
  gene_df <- data.frame(
    orig_name = names(expanded_gene_vector),
    Values = I(expanded_gene_vector)
  )
  # Expand the comma-separated values into multiple rows
  expanded_gene_df <- gene_df %>%
    separate_rows(Values, sep = ",\\s*") %>%
    rename(separated_gene_name = Values)
  
  return(expanded_gene_df)
}

get_DEGs <- function(sobj, ident="Genotype", cluster_name="post", df_gene_data) {
  # expect sobj to already be subsetted here to the group of interest
  #patient_num <- 6
  #sobj_PT_tumor <- subset(sobj, subset = PanCK.PT == 1 & Patient == patient_num & Run_Tissue_name %in% c("PSR-01", "PSR-02"))
  
  table(sobj$Sample.ID)
  
  Idents(sobj) <- ident
  sobj <- PrepSCTFindMarkers(sobj, assay = "SCT", verbose = TRUE)
  markers <- FindAllMarkers(sobj, assay="SCT", slot="data", only.pos = FALSE, verbose=TRUE, recorrect_umi = TRUE) # , min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  head(x = markers)
  print(paste("number of markers from FindAllMarkers", nrow(markers)))
  # Assuming you have already run FindAllMarkers() and have the markers dataframe
  # Add a column for significance threshold (e.g., adjusted p-value < 0.05 and log fold change > 0.25)
  markers <- markers %>%
    mutate(HighFoldChangeLowPval = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 1.0, "Y", "N")) %>%
    mutate(abs_logFC = abs(avg_log2FC))
  markers.pre <- markers %>% filter(cluster == cluster_name)
  head(markers.pre)
  
  table(markers.pre$HighFoldChangeLowPval)
  
  # clean up rownames to not have .1 at end. 
  rownames(markers.pre) <- gsub("\\.1", "", rownames(markers.pre))
  markers.pre$Gene <- rownames(markers.pre)
  
  # get gene annotations
  df_degs <- merge(markers.pre, df_gene_data, by.x = "Gene", by.y = "Display_Name", all.x = TRUE)
  df_degs <- df_degs %>% 
    dplyr::arrange(desc(abs_logFC))
  head(df_degs)
  # do this yourself.
  #write_csv(df_degs, paste0(results_dir, "/DEGs_PanCKTumor_regions_Patient",patient_num,".csv"))
  
  #df_degs <- df_degs %>% 
  #  dplyr::select(-pct.1, -pct.2, -gene, -cluster)

  return(df_degs)
}

plot_volcano_degs <- function(df_degs, num_labels = 50, title="", results_dir, 
                              plot_height=10, plot_width=10) {
  
  # Expect a column names "HighFoldChangeLowPval" with values "Y" or "N"
  # top_genes <- df_degs %>%
  #   dplyr::filter(HighFoldChangeLowPval == "Y") 
  top_genes <- df_degs %>%
    dplyr::top_n(num_labels, wt = abs(avg_log2FC))
  
  # Cap very low p-values at a threshold (e.g., 1e-300)
  df_degs$p_val_adj <- pmax(df_degs$p_val_adj, 1e-300)
  
  p <- ggplot(df_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = HighFoldChangeLowPval)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Y" = "red", "N" = "gray")) +
    #theme_minimal() +
    labs(title = paste("DEGs",title),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA)) +
    
    guides(color = "none") +   # Hide the color legend
    # Add labels for the top 10 highest fold change genes
    #coord_cartesian(ylim = c(0, 50)) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = Gene), 
                             size = 3, 
                             box.padding = 0.3, 
                             max.overlaps = Inf)
  
    print(p)
    ggplot2::ggsave(paste0(results_dir, "/VolcanoPlot_", sanitize_name(title), ".png"), 
                    plot=p, 
                    dpi = 300)
    return(p)
}

plot_volcano_degs_test <- function(df_degs, num_labels = 50, title="", results_dir, 
                              plot_height=10, plot_width=10) {
  
  # Expect a column names "HighFoldChangeLowPval" with values "Y" or "N"
  # top_genes <- df_degs %>%
  #   dplyr::filter(HighFoldChangeLowPval == "Y") 
  top_genes <- df_degs %>%
    dplyr::top_n(num_labels, wt = abs(avg_log2FC))
  
  # Cap very low p-values at a threshold (e.g., 1e-300)
  df_degs$p_val_adj <- pmax(df_degs$p_val_adj, 1e-300)
  
  p <- ggplot(df_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = HighFoldChangeLowPval)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Y" = "red", "N" = "gray")) +
    #theme_minimal() +
    labs(title = paste("DEGs",title),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA)) +
    
    guides(color = "none") +   # Hide the color legend

  print(p)
  return(p)
}



plot_deseq_volcano_degs <- function(df_degs, num_labels = 50, title="", results_dir) {
  
  top_genes <- df_degs %>%
    dplyr::top_n(num_labels, wt = abs(log2FoldChange))
  
  # Cap very low p-values at a threshold (e.g., 1e-300)
  df_degs$padj <- pmax(df_degs$padj, 1e-300)
  
  p <- ggplot(df_degs, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Y" = "red", "N" = "gray")) +
    theme_minimal() +
    labs(title = paste("Differentially Expressed Genes", title),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA)) +
    # Add labels for the top 10 highest fold change genes
    #coord_cartesian(ylim = c(0, 100)) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             size = 3, 
                             box.padding = 0.3, 
                             max.overlaps = Inf)
  
  print(p)
  ggplot2::ggsave(paste0(results_dir, "/DESeqVolcanoPlot_", sanitize_name(title), ".png"), plot=p,width = 8, height = 6, dpi = 300)
  return(p)
}

fetch_gene_data <- function(sobj, genes_list){
  # given a subset object and list of genes, fetch normalized data
  # and metadata needed for typical boxplots
  norm_data <- FetchData(sobj, gene_names, layer = "data")
  df_meta <- sobj@meta.data[c('TimePoint','Sample.ID','Sample.Label', 'PatientID')]
  gene_data <- merge(norm_data, df_meta, how="inner", by="row.names")
  return (gene_data)
  
}




plot_molecules <- function (sobj, fov_name, mol_list, mol_colors, title, show_legend=TRUE, molsize=0.1) {
  
  if (show_legend == TRUE) {
    legend.pos = "right"
  }
  else {
    legend.pos = "none"
  }
  p <- ImageDimPlot(sobj, fov = fov_name, axes = FALSE, # border.color = "white", border.size = 0.1, 
                    cols = "polychrome",
                    coord.fixed = TRUE, 
                    flip_xy = FALSE,
                    molecules = mol_list, 
                    nmols = 20000,
                    #crop=TRUE,
                    mols.cols = mol_colors, 
                    mols.size = molsize, mols.alpha = 0.5,
                    dark.background = FALSE) +
    ggtitle(title) +
    #scale_y_reverse() + 
    theme_minimal() +
    theme(legend.title = element_blank(),
          text = element_text(size = 10),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = legend.pos
          ) +
    coord_fixed(ratio=1)
  
  return (p)
  
}

niche_plot_spleen <- function(sobj, fov_name, colname, title  ) {
  
  Idents(sobj) <- colname
  p <- ImageDimPlot(sobj, fov = fov_name, axes = FALSE, # border.color = "white", border.size = 0.1, 
                    cols = "polychrome",
                    coord.fixed = TRUE, 
                    flip_xy = FALSE,
                    dark.background = FALSE) +
    ggtitle(title) +
    #scale_y_reverse() + 
    theme_minimal() +
    theme(legend.title = element_blank(),
          text = element_text(size = 10),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
    ) +
    coord_fixed(ratio=1)
  
  return(p)
}
