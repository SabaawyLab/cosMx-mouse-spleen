###########################
## ---------------------------
##
## Script name: mouse_spleen_RNA_EDA.R
##
## Purpose of script: load the Seurat data object produced by the AtoMx pipeline,
##   generate plots, evaluate rna quality controls
##.  summarize and display what has already been done in AtoMx
## 1. Use Seurat LoadNanostring() method (vs exported Seurat RDS) but either is fine. 
##    1b. subtract Negative Probes from the counts matrix
##    1c. Apply QC filters to the Seurat object
## 2. Load Sample metadata per cell:
##    2a. FOV and Sample.ID assignments
##    2b. Apply panCK/CD45 phenotypes to cells where possible. 
##    2c. overlay the pathology annotations to get additional meta: path_annot_id and name
## 3. Apply Cell Typing to the Seurat object with InSituType 
## 4. Split sobj into tissue, patient, and samples vs Organoids
## 5. Run PCA, UMAP on the Seurat object(s)
## 6. Cell Typing by Clustering
## 7. Create sub-fovs to visualize each Sample
## 8. plots and visualizations
## 10. Render array of plots to aid in selecting the most biologically correct cell typing/clustering
##
##
## Author: Ann Strange
##
## Date Created: 2024-10-14
##
## Email: ann.strange@cuanschutz.edu
##
## ---------------------------
##
## Notes:
##   What does our Seurat object from AtoMx look like? 
##.  especially the rna aspects. 
##. Key files:
##   
##
## ---------------------------
# 
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("huayc09/SeuratExtend")
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# install.packages("ggplot")
# install.packages("terra")
# devtools::install_version("Matrix", version = "1.6-5")
# if(!"pak" %in% installed.packages()) {
#   install.packages("pak")
# }
# install.packages('magick', repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))
# 
# #URL 'https://cran.rstudio.com/bin/macosx/big-sur-x86_64/contrib/4.4/spatstat.geom_3.3-2.tgz': Timeout of 60 seconds was reached
# #Warning in download.packages(pkgs, destdir = tmpd, available = available,  :
# #                               download of package ‘spatstat.geom’ failed
# install.packages("spatstat.geom")
# # stringi failed
# # RSpectra
# # shiny
# 
# pak::pkg_install("drieslab/Giotto")
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # Install Insitutype:
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("SummarizedExperiment")
# devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
# 
# #ERROR: dependencies ‘SingleCellExperiment’, ‘SummarizedExperiment’ are not available for package ‘InSituType’
# 
# 
# BiocManager::install("biomaRt")
# BiocManager::install("DESeq2")

#setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
#install.packages(c("BPCells", "presto", "glmGamPoi"))


require(tidyverse)
require(data.table)
library(Seurat)
library(Matrix)
library(Giotto)
#library(viridis)
#library(biomaRt) # gene alias lookups
library(ggplot2)
library(InSituType)
#library(trendsceek)
#library(multinet)
#library(RTriangle)
library(magick)
library(cowplot)
library(patchwork)
library(grid)
library(readxl)
library(pheatmap)
library(SeuratExtend)
library(DESeq2)
library(ggrepel)
library(openxlsx)
library(gridExtra)
#library(biomaRt)

old.wd <- getwd()           # save the current working directory
setwd("~/Documents/git/cosMx-mouse-spleen/")   # set working directory (mac)
source("R/mouse_seurat_utils.R") # load up the functions we need
#source("R/image_reg_tools.R")
#source("R/gene_annotations.R") # needs a little fixing
#source("R/cell_typing_utils.R")

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
#memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on macs.
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize= 10*1024*1024^2) # 10GB

## ---------------------------

## load up the packages we will need:  (uncomment as required)


# source("functions/packages.R")       # loads up all the packages we need

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------

###################################
# File structure setup, Mac syntax
###################################
# Modify these for each RNA experiment after export from AtoMx

#root_dir <- '/Volumes/T7Shield/PSR-GEN-057' # Experiment specific folder structure
root_dir <- '/Users/annmstrange/Documents/Projects/cosMx_MouseOct2024/Sabaawy_MuTMA_RNA_07_10_2024_21_25_12_667'
root_ref_dir <- '/Volumes/T7Shield' # reference dir for files not specific to this Experiment. 
exp_name <- 'Sabaawy_MuTMA_RNA_07_10_2024_21_25_12_667'
slide_names <- c( 'Sabaawynewcore091320246', 'Sabaawynewcore091320245')
# The cell ids in these are labeled c_1, c_2, c_3, c_4 respectively in the seurat object and metadata files.
slide_names
slide_name <- 'Sabaawynewcore091320245' # select one to work with (can be arbitrary)
# run_names is not needed but could be used to navigate the directory structure
run_names <- c('20240920_233038_S2','20240920_233038_S1')
run_name <- '20240920_233038_S1' # matches slide_name

# Caution: There is a cell naming convention in the flatfile vs the exported Seurat Object that differ. 
# After seurat::LoadNanostring(), we will rename the keys to match the Seurat object, like c_<slide_num>_<fov>_<cell_num> e.g. c_1_1_300
# cell_label_names is not explicitly needed, but merge will prepend these to the rownames if uniqueness is needed. 
# Our use of rename_keys makes the use of these cell_label_names unnecessary on merge.
# cell_label_names <- c('c_1', 'c_2', 'c_3', 'c_4')

results_dir <- paste0(root_dir, "/spleen_results")
# results_sub_wholeassay <- "/WholeAssay"
# results_sub_patient <- "/ByPatient"
# results_sub_region <- "/ByRegion"
# results_sub_celltype <- "/ByCellType"
# results_sub_panck <- "/ByPanCK"
# results_sub_annot <- "/ByAnnotatedRegion"
# results_sub_annot_celltype <- "/ByAnnotAndCellType"

rna_root_dir <- paste0(root_dir, '/RawFiles') # or flatFiles
rna_flatfiles_dir <- paste0(root_dir, '/flatFiles')
rna_dir <- paste0(rna_root_dir, '/', slide_name, '/', run_name, '/AnalysisResults/s2wi20fpqr') 
# This text file lists all the Probes
# /Users/annmstrange/Documents/Projects/cosMx_MouseOct2024/Sabaawy_MuTMA_RNA_07_10_2024_21_25_12_667/RawFiles/Sabaawynewcore091320245/20240920_233038_S1
rna_lookup <- paste0(rna_root_dir, '/', slide_name, '/', run_name, '/plex-s2wi20fpqr.txt')

metadata_dir <- paste0(root_dir,'/Metadata')
df_sample <- read.csv(paste0(metadata_dir, '/Sample_Metadata.csv'))
df_fov_meta <- read.csv(paste0(metadata_dir, '/FOV_Sample_Ids_RNA.csv'))

# This nanostring-provided file lists Probe names, gene names, and some annotations
markers_dir <- paste0(root_dir, '/BioMarkers')
ns_gene_data <- 'LBL-11176-03-Mouse-Universal-Cell-Characterization-Gene-List.xlsx'
ns_gene_fn <- paste0(markers_dir, '/', ns_gene_data)
df_gene_data <- get_nanostring_mouse_gene_annotations(ns_gene_fn )
head(df_gene_data, 3)

# Convert Human to Mouse
# load cosMx Cell Profile for Human IO
human_cell_profile <- read.csv( "/Users/annmstrange/Documents/git/CosMx-Cell-Profiles/Human/IO/IO.profiles.csv",
                                 row.names = 1)

head(human_cell_profile)

# inner join on Human_Gene
human_cell_profile$gene <- rownames(human_cell_profile)
new_matrix <- merge(human_cell_profile, df_gene_data, by.x ="gene", by.y= "Human_Gene")
head(new_matrix)
dim(human_cell_profile)
dim(new_matrix)
# Lost 100, which ones?
setdiff(human_cell_profile$gene, new_matrix$gene)
setdiff(df_gene_data$Display_Name, new_matrix$Display_Name)

human_genes <- rownames(human_cell_profile)
human_genes
new_matrix <- new_matrix[,1:17]
rownames(new_matrix) <- new_matrix$Display_Name
new_matrix <- new_matrix[, -1]
new_matrix <- new_matrix[, -16]
head(new_matrix)

df2 <- as.data.frame(lookup_list)

# get genes from linked databases
lookup_list <- convertHumanGeneList(c("TP53", "BRCA1", "EGFR"))
lookup_list <- convertHumanGeneList(human_genes)
length(lookup_list)

# given a list of Murine gene equivalents, there will be some gaps and dupes.
# now match each human gene to the mouse equivalent (NS Display Name) found for our Mouse experiment. 
head(ns_gene_data)
lookup_list2 <- findMatchingDisplayGene(c("Trp53", "Brca1"), ns_gene_data)


df_human_to_mouse <- data.frame(human = names(lookup_list), mouse = as.character(lookup_list))
#df_human_to_mouse <- as.data.frame(lookup_list)
write.csv(df_human_to_mouse, file = "/Users/annmstrange/Documents/git/CosMx-Cell-Profiles/Human/IO/human_to_mouse.csv")

# add to CellProfileMtx
human_cell_profile$mouse_genes <- mus_genes




# drop ns_annotations
#df_gene_data <- df_gene_data %>% dplyr::select(-ns_annotations)

df_special_genes <- read.csv(paste0(markers_dir, '/GenesOfInterest_RNA_HES.csv'))

colnames(df_special_genes) <- c('Display_Name', 'OfInterest_Reason', "Notes", "GroupAs", "Highlight")
df_special_genes <- df_special_genes %>% dplyr::select(Display_Name, OfInterest_Reason, Notes, GroupAs, Highlight)
head(df_special_genes)

genes_special <- df_special_genes$Display_Name
genes_special

genes_highlight <- df_special_genes %>%
  filter(Highlight == "x") %>%
  pull(Display_Name)

genes_highlight
# breakdown
genes_senescence <- df_special_genes %>%
  filter(OfInterest_Reason %in% c("Senescence")) %>%
  pull(Display_Name)

genes_special_il <- df_special_genes %>%
  filter(OfInterest_Reason %in% c("Interleukins", "Chemokines")) %>%
  pull(Display_Name)
genes_special_il

genes_special_gf <- df_special_genes %>%
  filter(OfInterest_Reason %in% c("Growth Factors")) %>%
  pull(Display_Name)

genes_bcells <- df_special_genes %>%
  filter(OfInterest_Reason %in% c("Bcell")) %>%
  pull(Display_Name)
genes_bcells 

# genes_special_ne <- c("S100B" )
# genes_p53 <- c("TP43", "RB1")
# 
# # don't bother plotting 
# genes_low_tumor_expr <- c("KRT5", "KRT13")

# left outer join with df_gene_data by RNA1K to gene 
df_gene_data2 <- left_join(df_gene_data, df_special_genes, by ='Display_Name')
df_gene_data2$OfInterest_Reason[is.na(df_gene_data2$OfInterest_Reason)] <- " "
df_gene_data2$Notes[is.na(df_gene_data2$Notes)] <- " "
df_gene_data2$GroupAs[is.na(df_gene_data2$GroupAs)] <- " "
df_gene_data2$Highlight[is.na(df_gene_data2$Highlight)] <- " "
head(df_gene_data2)
nrow(df_gene_data2)
df_gene_data <- df_gene_data2

table(df_gene_data2$OfInterest_Reason)

file.exists(results_dir)
file.exists(rna_root_dir) 
file.exists(rna_dir)
file.exists(rna_lookup)
file.exists(metadata_dir)
file.exists(paste0(metadata_dir, '/Sample_Metadata.csv'))
file.exists(paste0(metadata_dir, '/FOV_Sample_Ids_RNA.csv'))

#################################
# Resume here
#################################
saved_rds_dir <- paste0(root_dir, "/saved_objs")
file.exists(saved_rds_dir)

#sobj_load <- readRDS(paste0(saved_rds_dir, "/latest_sobj_15Oct24.rds"))
#sobj <- sobj_load


##################################
# Marker lists
#################################
#iocolors

# Leverage what ns did with iocolors for our final cell_types list, like this: 
iocolors <- c("B.cell"= "darkblue", 
                "Endothelial"=    "#996633" ,
                "Fibroblast" = "#999999" ,
                "Macrophage" =  "#006600" ,
                "Mast.cell" =  "springgreen",  # "#FFFF00" ,
                "mDC" = "#00FF00"  ,
                "Monocyte" = "#33CC00"  ,
                "Neutrophil" =  "#9966CC"  ,
                "NK.cell" =    "pink", # "grey20" ,
                 "pDC" = "#8DD3C7", #00FFFF" ,
                "DC" = "#8DD3C7", # "grey30" ,
                "Plasmablast" =   "#3399CC"  ,
                "Plasma-cell" = "#3399CC" ,
                "T CD4 memory" =  "skyblue", # "#FF0000"  ,
                "T CD4 naive" =  "turquoise", #CC0000",
                "T CD8 memory" =    "blue", #  "#FF9900" ,
                "T CD8 naive" =   "lightblue" , #   "#FF6633",
                "T-cell" = "#8DD3C7",
                "Treg" =   "cyan" 
              ) #  "#FF66FF" )
 
# meanwhile, unnamed color list, or polychrome
# celltype_colors <- c("darkblue",  "#996633" ,"#999999" ,"#006600" ,"#FFFF00", '#B3DE69', 
# "#33CC00"  ,"#9966CC"  ,"grey10", "#3399CC"  ,"#FF0000"  , "#CC0000", "lightblue",
#  "#FF9900" , "#FF6633", "#FF66FF" , 'purple', '#8DD3C7', '#FB8072', '#80B1D3',
# '#FDB462', 'magenta', 'blue', '#BC80BD','cyan', 'green' )

# celltype_colors <- c("B.cell"= "darkblue", 
#               
#               "Fibroblast" = "#999999" ,
#               "Macrophage" =  "#006600" ,
#               "Mast" = "#FFFF00" ,
#               "mDC" = "#00FF00"  ,
#               "Monocyte" = "#33CC00"  ,
#               "Neutrophil" =  '#B3DE69', 
#               "NK" =    '#80B1D3', 
#               "pDC" = "#00FFFF" ,
#               "Plasmablast" =   "#3399CC"  ,
#               "T.cell" = 'skyblue',
#               "T.CD4.memory" =   'purple',
#               "T.CD4.naive" =   '#BC80BD',
#               "T.CD8.memory" =     "#996633" ,
#               "T.CD8.naive" =      "pink",
#               "Treg" =    "#FF66FF",
#               "Endothelial"= '#FDB462',   
#               "Tumor" = "#FF0000",



# Load RNA labels for RNA marker attributes

df_rna_lookup <- read.table(rna_lookup, header = TRUE, sep = "\t")
head(df_rna_lookup) # ProbeID e.g CPROT01266, DisplayName e.g. 'CD4'
# example lookup; RNA has no ProbeID, but DisplayName
df_rna_lookup[df_rna_lookup$CodeClass != 'Endogenous', c('DisplayName','CodeClass')] 
# CodeClass can be SystemControl, Negative, or Endogenous. 

print(table(df_rna_lookup$CodeClass))
# expect counts similar to this:
#    Endogenous      Negative SystemControl 
# 950            10           197 

#######################################
# Build table of aliases, or other marker attributes
#######################################
# 
# ns_gene_ref <- '/Volumes/T7Shield/PSR-GEN-057/rna/PSR-GEN-057_ns_gene_reference.csv'
# df_ns_gene_ref <- read.csv(ns_gene_ref)
# head(df_ns_gene_ref)
# # which genes have Cell.Type == '+'
# cell_type_genes <- as.character(df_ns_gene_ref[df_ns_gene_ref$Cell.Type == '+', c('HUGO.Symbol')])
# cell_type_genes[1:5]

# Cell Typing per paper: ___________







###################
# seurat_dir <- '/Volumes/T7Shield/PSR-GEN-057/rna'
# #sobj_rds <- readRDS(paste0(seurat_dir, '/03513d30-eb51-44fd-92e9-fd2ea8ec1890_seuratObject.RDS'))
# 

# Alternative Method to load Seurat object and metadata from RDS file: 
# The RDS Seurat object does offer the x and y coords in mm, but otherwise the metadata is very similar
# We can consider them nearly equivalent. 
#sobj_rds <- readRDS(paste0(seurat_dir, '/03513d30-eb51-44fd-92e9-fd2ea8ec1890_seuratObject.RDS'))


########################################
# Bugfix, invert all the tx coords orientations
#######################################

# If inverted, CenterY_global_px will be < 10K, otherwise > 125K
#df_metadata[df_metadata$cell_id == "c_3_1_1", c('CenterY_global_px', 'CenterX_global_px', 'fov','Run_Tissue_name')]


#########################################
# Seurat LoadMetadata() from flat files
#########################################

  
obj_list <- list()
negmeans_list <- list()

# nanostring case:
#flatfiles_dir <- paste0(rna_root_dir, '/', slide_names[1], '/', run_names[1])
#obj1_ex <- LoadNanostring(data.dir = flatfiles_dir, fov = slide_names[1], assay="Nanostring")
#saveRDS(obj1_ex, file = paste0(results_dir, '/obj1_postLoadNS.RDS'))

load_mouse_meta_load_missed <- function(sobj, metadata_file) {
  # Seurat LoadNanostring() workaround to load metadata flat file including making sure the keys match
  # this loads morph markers, measurements, X, Y, (global and local) and cell_ID
  # TODO: add error checking on cell_ID format as integer
  
  #df_obj6_meta <- read.csv(paste0(lung6.data.dir, "/Lung6_metadata_file.csv"))
  df_obj_meta <- read.csv(metadata_file)
  print(colnames(df_obj_meta))
  # workaround to fix inverted centroids, flips CenterY_global_px
  #df_obj_meta <- flip_centroids(df_obj_meta)
  
  # shorten unweildy column names
  string_to_replace <- "a57c6ac1.3d04.4093.af08.340b8f4c3e0d"
  replacement_string <- "a57"
  names(df_obj_meta) <- gsub(string_to_replace, replacement_string, names(df_obj_meta))
  string_to_replace <- "5f7f57b1.56ee.4536.b017.b5ec1fce7297"
  replacement_string <- "af7"
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

i=2
for (i in 1:2){
  print(i)
  print(slide_names[i])
  print(run_names[i])
  flatfiles_dir <- paste0(rna_flatfiles_dir, '/', slide_names[i])
  print(flatfiles_dir)
  print(file.exists(flatfiles_dir))
  list.files(flatfiles_dir)
  
  # Fix polygons file to re-orient global y coords, do only once. 
  #success <- reorient_polygon_global_y(flatfiles_dir, slide_name = slide_names[i])
  
  obj_list[[i]] <- LoadNanostring(data.dir = flatfiles_dir, fov = slide_names[i], assay="RNA")
  # LoadNanostring doesn't handle the metadata 
  obj_list[[i]] <- load_mouse_meta_load_missed(obj_list[[i]], paste0(flatfiles_dir, "/",slide_names[i],"_metadata_file.csv.gz"))
  # This rename is needed as LoadNanostring uses <cell>_<fov> cell ID vs the c_<slide>_<fov>_<cell> format used elsewhere
  obj_list[[i]]  <- rename_keys(obj_list[[i]]) 
  # Each slide should be normalized independently 
  
  # Remove SystemControls
  print(DefaultAssay(obj_list[[i]]))
  #Error in (function (cl, name, valueClass)  : 
  #            ‘counts’ is not a slot in class “Assay5”
  print(paste0("Num Features before removing System Controls: ", nrow(obj_list[[i]])))
  # this did nothing:
  #SetAssayData(obj_list[[i]], layer = "counts", new.data = remove_sys_control(obj_list[[i]], "Nanostring", "counts"))
  obj_list[[i]] = remove_sys_control(obj_list[[i]], "RNA", "counts")
  print(paste0("Num Features after removing System Controls (expect 1010): ", nrow(obj_list[[i]])))
  
  #obj_list[[i]][["Nanostring"]]@counts <- remove_sys_control(obj_list[[i]], "Nanostring", "counts")
  
  # Negative controls should be subtracted/removed before normalization, but not before running InSituType
  # However, we keep 0 values for Negative controls in the counts matrix after subtracting, and set aside the means for later
  negmeans_list[[i]] <- get_neg_control_means(obj_list[[i]], "RNA", "counts")
  # add metadata for neg
  obj_list[[i]] <- AddMetaData(obj_list[[i]], metadata = data.frame(neg = negmeans_list[[i]]))
 
  #SetAssayData(obj_list[[i]], layer = "counts", new.data = subtract_neg_control(obj_list[[i]], "Nanostring", "counts"))
  obj_list[[i]] <- subtract_neg_control(obj_list[[i]], "RNA", "counts")
  print(paste0("Num Features after accounting for Negative Controls (s/b unchanged): ", nrow(obj_list[[i]])))
  print(paste0("Num Cells before QC: ", ncol(obj_list[[i]])))
  # BASIC QC filter needed to prevent sparsity errors in SCTransform; option for more later 
  obj_list[[i]] <- subset(obj_list[[i]], subset = (nFeature_RNA > 10 & nCount_RNA > 20)) #20 genes/cell. also vs nFeature_Nanostring
  print(paste0("Num Cells after super basic QC: ", ncol(obj_list[[i]])))
  # Remove cells that didn't pass the AtoMx QC
  obj_list[[i]] <- subset(obj_list[[i]], subset = qcCellsFlagged == FALSE)
  print(paste0("Num Cells after filtering on more complete cell QC: ", ncol(obj_list[[i]])))
  
  # Normalizes counts can be done later after QC
  #obj_list[[i]] <- SCTransform(obj_list[[i]], assay = "Nanostring", verbose = TRUE)
  # This should match the Seurat method in AtoMx
  obj_list[[i]] <- NormalizeData(obj_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  obj_list[[i]]
  
}

Layers(obj_list[[i]])

# Merge Seurat objects using chaining
# Merging allows for group metadata to be added to the Seurat object, and also allows for
# joint dimensional redux and clustering on the underlying RNA expression data if desired. (but we won't)
sobj <- merge(obj_list[[1]], y = c(obj_list[[2]]))
#sobj <- obj_list[[i]]
head(sobj@meta.data$neg)
sum(sobj@meta.data$neg)

colnames(sobj@meta.data)

#class(negmeans_list[[1]])
#negmeans_vec <- c(negmeans_list[[1]], negmeans_list[[2]], negmeans_list[[3]], negmeans_list[[4]])
#sum(negmeans_vec) # 4002 
# Save for later
#write.table(negmeans_vec, file = paste0(results_dir, "/negmeans_list.txt"), row.names=TRUE, col.names = FALSE)

# Save the Seurat object 
saveRDS(sobj, file = paste0(saved_rds_dir, "/sobj_from_Load_v1.RDS"))
# to resume here after LoadNanostring
#sobj <- readRDS(paste0(saved_rds_dir, "/sobj_from_Load_v1.RDS"))


#################################################
# Add Sample and FOV Metadata
################################################


# set rownames
rownames(df_sample) <- df_sample$Sample.ID

head(df_sample)
# keep only subset of columns
df_sample <- df_sample[, c('Sample.ID', 'Sample.Label', 'Sample.Nm', 'tissue', 'condition', 'Physical.Tag', 'Sex', 'Date.of.Birth', 'Genotype', 'time.to.form.tumor', 'Organ')]

# calculate nicer looking labels
# df_sample <- df_sample %>%
#   mutate(Sample.Label = paste(paste(PatientID, tissue, TimePoint, sep = "-"),"(", Sample.ID, ")")) %>%
#   mutate(Sample.Label2 = paste0(paste(SPORE, tissue, TimePoint, sep = "-")," (", Sample.ID, ")")) 

head(df_fov_meta)
#df_fov_meta <- df_sample[, c('Run_Tissue_name', 'fov', 'Sample.ID')]
  

# Replace "- " with " " in the Sample.Label column
df_sample$Sample.Label <- gsub("\n", "", df_sample$Sample.Label)
df_sample$Sample.Label <- gsub("  ", " ", df_sample$Sample.Label)
#df_sample$Sample.Label2 <- gsub("- ", " ", df_sample$Sample.Label2)


#join 
df_fov_meta2 <- merge(df_fov_meta, df_sample, by = 'Sample.ID')
head(df_fov_meta2)
table(df_fov_meta2$Sample.Label)
table(df_fov_meta2$Genotype)

# df_fov_meta2 <- df_fov_meta2 %>%
#   mutate(Sample.Label3 = case_when(Diagnosis != "" ~ paste(Diagnosis, "cell line"), 
#                                    TRUE ~ Sample.Label2))
# unique(df_fov_meta2$Sample.Label3)
# df_fov_meta2$Sample.Label3 <- factor(df_fov_meta2$Sample.Label3, 
#                                     levels = c("220764-lung-post (1_1)",
#                                                "CUTO-43-lung-post (2_4)",
#                                                "CUTO-59-lymph node-pre (3_5)" ,
#                                                "CUTO-59-liver-post (3_6)",
#                                                "220904-liver-pre (4_6)",
#                                                "220904-liver-post (4_1)",
#                                                "220900-pleura-Rx-post (5_2)" ,
#                                                "220901-brain-pre (6_3)",
#                                                "220901-lymph node-post (6_5)",
#                                                "220899-T2-vertebrae-pre (7_4)",
#                                                "CUTO-43-PDO-post (ORG_2_3)",
#                                                "220904-PDO-post (ORG_4_1)", 
#                                                "CUTO-41-PDO-post (ORG_8_2)" ,
#                                                "220764-PDO-post (ORG_8_5)" ,
#                                                "Adenocarcinoma cell line",
#                                                "AdenoSquamous cell line",
#                                                "SCC cell line"))


#table(df_fov_meta2$Sample.Label3)

#sobj_tumor_pairs <- AddMetaData(df_fov_meta2, df_meta["Sample.Label3"])

# assign metadata to each cell by fov and flow cell (slide) i.e Run_Tissue_name
#df_metadata <- sobj@meta.data




sobj <- add_mouse_sample_metadata(sobj, df_fov_meta2)
table(sobj@meta.data$Organ)
table(sobj@meta.data$Genotype)
sobj

# Exclude FOVs marked for exclusion in the fov metadata
#table(sobj@meta.data$Exclude) # 2823 cells to remove
#sobj <- subset(sobj, subset = Exclude != "x")


# table(sobj@meta.data$Sample.ID, sobj@meta.data$PanCK.PT)
# table(sobj@meta.data$Sample.ID, sobj@meta.data$CD45.PT)

#sobj@meta.data$PanCK.PT[is.na(sobj@meta.data$PanCK.PT)] <- 1
#sobj@meta.data$CD45.PT[is.na(sobj@meta.data$CD45.PT)] <- 0

table(sobj@meta.data$Run_Tissue_name)
# what are the fovs?
sobj@images[1] # Sabaawynewcore091320246
Idents(sobj) <- "Sample.ID"
Idents(sobj) <- "condition"
# Todo: use subset of cells for performance
ImageDimPlot(sobj, fov="Sabaawynewcore091320246", col="glasbey")

Idents(sobj) <- "tissue"
ImageDimPlot(sobj, fov="Sabaawynewcore091320245")

Idents(sobj) <- "CD45.PT"
ImageDimPlot(sobj, fov="Sabaawynewcore091320246")
ImageDimPlot(sobj, fov="Sabaawynewcore091320245")

sobj_spleen <- subset(sobj, subset = Organ == "Spleen" & Sample.ID >= "30")
Idents(sobj_spleen) <- "Genotype"
# Todo: use subset of cells for performance
ImageDimPlot(sobj_spleen, fov="Sabaawynewcore091320246", col=c("red","blue","green"))


###################################
# Split things up
##################################
sobj_spleen <- subset(sobj, subset = Organ == "Spleen" & Sample.ID >= "30")
sobj_spleen5 <- subset(sobj_spleen, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen6 <- subset(sobj_spleen, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 6")



sobj_spleen_even <- subset(sobj, subset = Sample.ID %in% c(30,32,34)) 
sobj_spleen_even5 <- subset(sobj_spleen_even, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen_even6 <- subset(sobj_spleen_even, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 6")
                           
# rename to P14^WAF1/Cip1
feat_names <- rownames(sobj_spleen_even5)
feat_names <- gsub("Cdkn1a", "P21", feat_names)
rownames(sobj_spleen_even5) <- feat_names

table(sobj_spleen@meta.data$Sample.ID)
# Focusing on only 3 samples, slide 5
sobj_spleen30 <- subset(sobj_spleen5, subset = Sample.ID == 30 & 
                          Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen31 <- subset(sobj_spleen5, subset = Sample.ID == 31 & 
                          Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen32 <- subset(sobj_spleen5, subset = Sample.ID == 32 &
                          Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen34 <- subset(sobj_spleen5, subset = Sample.ID == 34 &
                          Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_spleen35 <- subset(sobj_spleen5, subset = Sample.ID == 35 & 
                          Run_Tissue_name == "Sabaawy new core 09/13/2024 5")

genes_senescence
#genes_senescence <- gsub("Cdkn1a", "P14", genes_senescence)

rownames(sobj_spleen30)[140:200]


features_genes <- c("P21", "Mki67" , "Lmna", "Chek1","Stmn1", "Pten", "Ccl5", "Cxcl13", "Gzmb" )
Idents(sobj_spleen5) <- "Genotype"
ImageDimPlot(sobj_spleen5, fov="Sabaawynewcore091320245")
# FeaturePlot
p <- ImageFeaturePlot(sobj_spleen30, fov="Sabaawynewcore091320245",
                 features = features_genes, 
                 size = 1,
                 coord.fixed=TRUE,
                 dark.background = FALSE) + 
  plot_layout(ncol=3)
p
# returns patchworked object

str(p)
ggsave(paste0(results_dir,"/spleen30_senescence_FeaturePlots.png"), plot = p,  width = 10, height = 12, dpi = 300)

p <- ImageFeaturePlot(sobj_spleen32, fov="Sabaawynewcore091320245",
                 features = features_genes, 
                 size = 1,
                 dark.background = FALSE) + 
    plot_layout(ncol=3)
ggsave(paste0(results_dir,"/spleen32_senescence_FeaturePlots.png"), plot = p, width = 10, height = 12, dpi = 300)


p <- ImageFeaturePlot(sobj_spleen34, fov="Sabaawynewcore091320245",
                 features = features_genes, 
                 size = 1,
                 dark.background = FALSE) + 
  plot_layout(ncol=3)
ggsave(paste0(results_dir,"/spleen34_senescence_FeaturePlots.png"), plot = p, width = 10, height = 12, dpi = 300)


print( theme_get())
#############################
# Sample Level Stats
#############################

# # TODO: move this below path annotations
# df_sample_stats <- sobj@meta.data %>%
#   group_by(Sample.ID) %>%  # Group by Sample.ID
#   summarise(
#     cell_count = n(),               # Count the number of rows per Sample.ID
#     unique_fovs = n_distinct(fov),   # Count the number of unique fovs
#     PanCK_tumor_count = sum(PanCK.PT == 1),
#     Path_annot_tumor_count = sum(path_annot_name.1 == "Annotated Tumor"),
#     #PanCK_pct_tumor = (PanCK_tumor_count / cell_count) * 100  # Percentage of tumor cells (PanCK.PT == 1)
#   )
# head(df_sample_stats)


#sobj_spleen <- subset(sobj, subset = Organ == "Spleen" & Sample.ID >= "30")
#table(sobj_spleen@meta.data$Sample.ID)
# Do we have any SCTransforms?
#obj_spleen@assays$RNA@SCTModel.list 

# normalize for spleen only
#sobj_spleen <- SCTransform(sobj_spleen, assay = "Nanostring", verbose = TRUE)

saveRDS(sobj_spleen_even, paste0(saved_rds_dir,"/sobj_spleen_even.rds"))
# get counts of each molecule by sample.

get_counts_subset_df <- function(sobj){
  
  exp_data <- as.matrix(GetAssayData(sobj, assay="RNA", layer="counts"))
  class(exp_data)
  dim(exp_data) # 1000 x 191K
  # how to subset matrix by gene
  
  return(as.data.frame(t(exp_data)))
}


exp_data <- as.matrix(GetAssayData(sobj_spleen_even, assay="RNA", layer="counts"))
norm_exp_data <- as.matrix(GetAssayData(sobj_spleen_even, assay="RNA", layer="data"))
class(exp_data)
dim(exp_data) # 1000 x 191K
# how to subset matrix by gene

# which genes in genes_special not found as rownames in exp_data
genes_special_not_found <- setdiff(genes_special, rownames(exp_data))
print(genes_special_not_found)
genes_special <- genes_special[!genes_special %in% genes_special_not_found]

exp_data2 <- exp_data[genes_special, ]
norm_exp_data2 <- norm_exp_data[genes_special,]
dim(exp_data2) # 94 x 191K

# percent of values in exp_data2 for each gene that are non-zero


percent_non_zero <- colMeans(exp_data2 > 0) * 100


#### plot lognormal data of n+1
log_data <- log(norm_exp_data2["Mki67", ]+1)
head(log_data)
hist(log_data, breaks = 100, col = "lightblue", 
main = "Mki67 Expression", xlab = "log10 Expression")

#########################
# Geometric Mean
# first show how to reverse this. 
original_data <- exp(log_data) - 1
head(original_data)

geometric_mean <- exp(mean(log_data, na.rm = TRUE))
geometric_mean

# Calculate the standard deviation on the log-transformed data
# After log transformation, the standard deviation reflects the variability 
# in the log-transformed scale rather than the original scale. If needed, you 
# can back-transform the results for interpretation on the original scale.
sd_log_transformed <- sd(log_data, na.rm = TRUE)
sd_log_transformed
sd_orig_scale <- exp(sd_log_transformed) 
sd_orig_scale

df_meta <- sobj_spleen_even@meta.data
# get cell counts per sample
cell_counts_by_sample <- table(df_meta$Sample.ID)

pct_non_zero <- function(x) {
  length(x[x >0]) / length(x) * 100
}

mean_non_zero <- function(x) {
  mean(x[x >0])
}


############################
# Cell Types
# handle niches separately  "spatialclust_af7_1_assignments"
############################
df_cell_types <- sobj_spleen_even@meta.data %>%
   dplyr::select ("cell", "Sample.ID", 
                  "RNA_nbclust_e3def5ab.413a.4c02.b683.f76ba4202d63_1_clusters",
                  #"RNA_nbclust_853a9cfe.66b9.4d5b.895c.daee922bc6a6_1_clusters",
                  "spatialclust_4786da16.29dc.44c6.ad08.97cef88cc2ee_1_assignments"
            )
  # group_by(Sample.ID, RNA_nbclust_853a9cfe.66b9.4d5b.895c.daee922bc6a6_1_clusters
  #         ) %>%
  # summarise(count = n(), .groups = "drop")
head(df_cell_types)
colnames(df_cell_types) <- c("cell", "Sample.ID", "CellType", "Niche")

# Assign "a" to "Macrophges" 
df_cell_types$CellType <- ifelse(df_cell_types$CellType == "c", "Macrophage", df_cell_types$CellType)

table(df_cell_types$CellType)
head(df_cell_types)


############################
# Check for Mutual Exclusivity
#############################

#df_meta <- sobj_spleen_even@meta.data
sobj_spleen30
# count matrix need
df_30 <- get_counts_subset_df(sobj_spleen30)
df_31 <- get_counts_subset_df(sobj_spleen31)
df_32 <- get_counts_subset_df(sobj_spleen32)
df_34 <- get_counts_subset_df(sobj_spleen34)
df_35 <- get_counts_subset_df(sobj_spleen35)

dim(df_30)
dim(df_31)
dim(df_32)

# Add phenotypes
df_30$cell <- rownames(df_30)
df_31$cell <- rownames(df_31)
df_32$cell <- rownames(df_32)
df_34$cell <- rownames(df_34)
#df_30$prolif <- ifelse(df_30$Pcna > 0 | df_30$Mki67 > 0, 1, 0)
#df_30$sen <- ifelse(df_30$Cdkn1b > 0 , 1, 0)



df_30 <- df_30 %>%
  mutate(prolif = ifelse(Pcna > 0 | Mki67 > 0, 1, 0)) %>%
  mutate(sen = ifelse(Cdkn1a > 0, 1, 0)) %>%
  mutate(sen_pt = ifelse(prolif == 1 & sen ==1, "both",
                  ifelse (prolif == 1, "prolif",
                    ifelse(sen == 1, "sen", "neither"))))
                    
table(df_30$sen_pt)

df_32 <- df_32 %>%
  mutate(prolif = ifelse(Pcna > 0 | Mki67 > 0, 1, 0)) %>%
  mutate(sen = ifelse(Cdkn1a > 0, 1, 0)) %>%
  mutate(sen_pt = ifelse(prolif == 1 & sen ==1, "both",
                         ifelse (prolif == 1, "prolif",
                                 ifelse(sen == 1, "sen", "neither"))))
table(df_32$sen_pt)

df_34 <- df_34 %>%
  mutate(prolif = ifelse(Pcna > 0 | Mki67 > 0, 1, 0)) %>%
  mutate(sen = ifelse(Cdkn1a > 0, 1, 0)) %>%
  mutate(sen_pt = ifelse(prolif == 1 & sen ==1, "both",
                         ifelse (prolif == 1, "prolif",
                                 ifelse(sen == 1, "sen", "neither"))))
table(df_34$sen_pt)
# add cell type to these, lookup for df_cell_types
df_30 <- merge(df_30, df_cell_types, by="cell")
dim(df_30)
table(df_30$CellType)

# celltype breakdown for each
df_30_summary <- df_30 %>%
  group_by(sen_pt, CellType) %>%
  summarize(count = n(), .groups = "drop") %>%
  #group_by(sen_pt) %>%
  group_by(CellType) %>%
  mutate(pct_celltype_inPT = (count / sum(count)) * 100) 

# df_30_summary <- df_30_summary %>%
#   dplyr::filter(percentage > 5) %>%
#   dplyr::filter(sen_pt != "neither")
df_30_summary$Sample.ID <- 30

df_32 <- merge(df_32, df_cell_types, by="cell")
dim(df_32)
# celltype breakdown for each
df_32_summary <- df_32 %>%
  group_by(sen_pt, CellType) %>%
  summarize(count = n(), .groups = "drop") %>%
  #group_by(sen_pt) %>%
  group_by(CellType) %>%
  mutate(pct_celltype_inPT = (count / sum(count)) * 100) 
# pct of the cell type in the PT (e.g. sen or prolif or both or neither)

# df_32_summary <- df_32_summary %>%
#   dplyr::filter(percentage > 5) %>%
#   dplyr::filter(sen_pt != "neither")
df_32_summary$Sample.ID <- 32

df_34 <- merge(df_34, df_cell_types, by="cell")
dim(df_34)
# celltype breakdown for each
df_34_summary <- df_34 %>%
  group_by(sen_pt, CellType) %>%
  summarize(count = n(), .groups = "drop") %>%
  #group_by(sen_pt) %>%
  group_by(CellType) %>%
  mutate(pct_celltype_inPT = (count / sum(count)) * 100) 

# df_34_summary <- df_34_summary %>%
#   dplyr::filter(percentage > 5) %>%
#   dplyr::filter(sen_pt != "neither")
df_34_summary$Sample.ID <- 34

# rbind
df_spleen_sen_summary <- rbind(df_30_summary, df_32_summary, df_34_summary)
df_spleen_sen_summary <- df_spleen_sen_summary %>%
  filter(count > 20) %>%
  filter(CellType %in% c( "B.cell", "c", "Macrophage", "Monocyte","Plasma" ,         
                          "Plasmablast" ,  "T.cell.CD4" , "T.cell.CD8","T.cell.regulatory"))

# pie plot for each sample
df_spleen_sen_summary$sen_pt <- factor(df_spleen_sen_summary$sen_pt, levels = c("both", "prolif", "sen", "neither"))
df_spleen_sen_summary$CellType <- factor(df_spleen_sen_summary$CellType)
df_spleen_sen_summary <- df_spleen_sen_summary %>%
  filter(sen_pt != "neither")

# levels = c("B cell", "T cell", "Macrophage", "DC", "Neutrophil", "Monocyte", "NK cell", "Mast cell", "Erythrocyte", "Plasma cell", "Basophil", "Megakaryocyte
# "))
# ggplot(df_spleen_sen_summary, aes(x="", y=percentage, fill=CellType)) +
#   geom_bar(width = 1, stat = "identity") +
#   coord_polar("y", start=0) +
#   facet_wrap(~Sample.ID + sen_pt) +
#   theme_minimal() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         legend.position = "bottom") +
#   labs(title = "Spleen Cell Type Breakdown by sen_pt",
#        fill = "Cell Type",
#        y = "Percentage",
#        x = "Sample ID")

levels(df_spleen_sen_summary$CellType)

# Barplot grouped by sample ID and sen_pt
ggplot(df_spleen_sen_summary, aes(x=sen_pt, y=pct_celltype_inPT, fill=CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Sample.ID) +
  theme_classic() +
  theme( #xis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  labs(title = "Spleen Cell Type Breakdown by Senescent/Prolif Phenotypes",
       fill = "CellType",
       y = "percent of CellType",
       x = "Phenotype")
ggsave(paste0(results_dir, "/spleen_sen_pt_celltype_barplot.png"), width = 10, height = 6, dpi=300)

# Filter columns based on the condition
num_cells_both30 <- df_30 %>%
  filter((Pcna > 0 | Mki67 > 0) & Cdkn1a > 0) %>%
  nrow()
num_cells_both30

num_cells_prolif30 <- df_30 %>%
  filter(Pcna > 0 | Mki67 > 0) %>%
  nrow()
num_cells_prolif30
num_cells_sensec30 <- df_30 %>%
  filter(Cdkn1a > 0) %>%
  nrow()
num_cells_sensec30
tot_cells30 <- nrow(df_30)

num_cells_both32 <- df_32 %>%
  filter((Pcna > 0 | Mki67 > 0) & Cdkn1a > 0) %>%
  nrow()
num_cells_both32

num_cells_prolif32 <- df_32 %>%
  filter(Pcna > 0 | Mki67 > 0) %>%
  nrow()
num_cells_prolif32
num_cells_sensec32 <- df_32 %>%
  filter(Cdkn1a > 0) %>%
  nrow()
num_cells_sensec32
tot_cells32 <- nrow(df_32)

num_cells_both34 <- df_34 %>%
  filter((Pcna > 0 | Mki67 > 0) & Cdkn1a > 0) %>%
  nrow()
num_cells_both34

num_cells_prolif34 <- df_34 %>%
  filter(Pcna > 0 | Mki67 > 0) %>%
  nrow()
num_cells_prolif34
num_cells_sensec34 <- df_34 %>%
  filter(Cdkn1a > 0) %>%
  nrow()
num_cells_sensec34
tot_cells34 <- nrow(df_34)
# need total cell counts by sample


# put in a dataframe 
sample.ids <- c(30,32,34)
both_counts <- c(num_cells_both30, num_cells_both32, num_cells_both34)
prolif_counts <- c(num_cells_prolif30, num_cells_prolif32, num_cells_prolif34)
sene_counts <- c(num_cells_sensec30, num_cells_sensec32, num_cells_sensec34)
pt_df <- data.frame(sample.id = sample.ids, 
                    both = both_counts, 
                    prolif = prolif_counts, 
                    sene = sene_counts)
pt_df

gene <- "Mki67"
# Step 4: Calculate mean and standard deviation per sample for each gene
gene_stats <- lapply(genes_special, function(gene) {
  # Extract expression for the gene
  gene_expression <- exp_data2[gene, ]
  norm_expression <- norm_exp_data2[gene, ]
  # Combine expression with sample information
  gene_df <- data.frame(expression = gene_expression,
                        norm_expression = norm_expression,
                        sample = df_meta$Sample.ID,
                        label = df_meta$Sample.Label,
                        genotype = df_meta$Genotype)
  
  # Calculate mean and standard deviation by sample
  # stats <- aggregate(expression ~ sample, data = gene_df, 
  #                    FUN = function(x) c(mean = mean(x), sd = sd(x)))
  # norm_stats <- aggregate(norm_expression ~ sample, data = gene_df, 
  #                    FUN = function(x) c(mean = mean(x), sd = sd(x)))
  
  # by sample
  stats <- aggregate(expression ~ sample, data = gene_df, 
                     FUN = function(x) c(pct_non_zero = pct_non_zero(x), mean_non_zero = mean_non_zero(x)))
  
  norm_stats <- aggregate(norm_expression ~ sample, data = gene_df, 
                      FUN = function(x) c(pct_non_zero = pct_non_zero(x), mean_non_zero = mean_non_zero(x)))
  
  stats
  norm_stats
  # non_zero_counts <- gene_expression[gene_expression >0]
  # print(length(non_zero_counts))
  # non_zero_counts_norm <- norm_expression[norm_expression >0]
  # print(length(non_zero_counts_norm))  # expected to be the same
  # pct_non_zero <- length(gene_expression[gene_expression >0]) / length(gene_expression) * 100
  # mean_pos_expr <- mean(non_zero_counts)
  # mean_pos_expr_norm <- mean(non_zero_counts_norm)
  
  # Add the number of cells for each sample to the result
  stats$cell_count <- cell_counts_by_sample[stats$sample]
  
  # join stats and norm_stats
  stats <- merge(stats, norm_stats, by="sample")
  #stats <- merge(pct_non_zero, mean_pos_expr, mean_pos_expr_norm, by="sample")
  # Format the result
  # stats <- data.frame(sample = stats$sample,
  #                     mean = sapply(stats$expression, "[", 1),
  #                     sd = sapply(stats$expression, "[", 2))
  stats$gene <- gene
  stats
})
df_gene_stats <- do.call(rbind, gene_stats)
class(df_gene_stats)
df_gene_stats
colnames(df_gene_stats) #<- c("sample", "expression.mean", "expression.sd", "gene"))
df_gene_stats <- df_gene_stats[, c("gene", "sample","cell_count", "expression", "norm_expression")]
head(df_gene_stats)

###############################
# GROUP BYs
###############################
# get unique GroupAs from df_gene_data2
grouped_genes <- unique(df_gene_data2$GroupAs) # CD3, CD8, VEGF
grouped_genes <- grouped_genes[!grouped_genes %in% c("", " ")]
grouped_genes
group <- "CD3"

# get stats by GroupAs
gene_stats_grouped <- lapply(grouped_genes, function(group) {
  # Extract expression for the gene
  gene_expression <- colSums(exp_data2[df_gene_data2[df_gene_data2$GroupAs == group, "Display_Name"], ])
  norm_expression <- colSums(norm_exp_data2[df_gene_data2[df_gene_data2$GroupAs == group, "Display_Name"], ])
  # Combine expression with sample information
  # not used
  gene_df <- data.frame(expression = gene_expression,
                        norm_expression = norm_expression,
                        sample = df_meta$Sample.ID,
                        label = df_meta$Sample.Label,
                        genotype = df_meta$Genotype
                        )
  
  # Calculate mean and standard deviation by sample
  # stats <- aggregate(expression ~ sample, data = gene_df, 
  #                    FUN = function(x) c(mean = mean(x), sd = sd(x)))
  # norm_stats <- aggregate(norm_expression ~ sample, data = gene_df, 
  #                    FUN = function(x) c(mean = mean(x), sd = sd(x)))
  # by sample
  stats <- aggregate(expression ~ sample, data = gene_df, 
                     FUN = function(x) c(pct_non_zero = pct_non_zero(x), mean_non_zero = mean_non_zero(x)))
  
  norm_stats <- aggregate(norm_expression ~ sample, data = gene_df, 
                          FUN = function(x) c(pct_non_zero = pct_non_zero(x), mean_non_zero = mean_non_zero(x)))
  
  # join stats and norm_stats
  stats <- merge(stats, norm_stats, by="sample")
  
  # Add the number of cells for each sample to the result
  stats$cell_count <- cell_counts_by_sample[stats$sample]
  
  # Format the result
  # stats <- data.frame(sample = stats$sample, 
  #                     mean = sapply(stats$expression, "[", 1), 
  #                     sd = sapply(stats$expression, "[", 2))
  stats$gene <- group
  stats
})

gene_stats_grouped <- do.call(rbind, gene_stats_grouped)
class(gene_stats_grouped)
gene_stats_grouped
colnames(gene_stats_grouped) #<- c("sample", "expression.mean", "expression.sd", "gene"))
gene_stats_grouped <- gene_stats_grouped[, c("gene", "sample","cell_count", "expression", "norm_expression")]
head(gene_stats_grouped)

# Now combine the two 
colnames(gene_stats_grouped) #<- c("sample", "expression.mean", "expression.sd", "gene"))
gene_stats_grouped <- gene_stats_grouped[, c("gene", "sample","cell_count", "expression", "norm_expression")]
colnames(df_gene_stats)
df_list <- list(df_gene_stats, gene_stats_grouped)
df_gene_stats2 <- do.call(rbind, df_list)
 

# add gene names from df_gene_data2
df_gene_data_temp <- df_gene_data2[,c("Display_Name", "OfInterest_Reason", "Notes", "GroupAs", "Highlight")]
df_gene_stats3 <- merge(df_gene_stats2, df_gene_data_temp, by.x="gene", by.y="Display_Name")
head(df_gene_stats3)

# also get sample.label from df_sample
df_gene_stats3 <- merge(df_gene_stats3, df_sample[,c("Sample.ID", "Sample.Nm", "Genotype")], by.x="sample", by.y="Sample.ID")
head(df_gene_stats3)

# Add calculations for SEM to each column
#df_gene_stats3$exp_stderr_ofmean  <- df_gene_stats3$expression[,"sd"] / sqrt(df_gene_stats3$cell_count)
#df_gene_stats3$norm_stderr_ofmean <- df_gene_stats3$norm_expression[,"sd"] / sqrt(df_gene_stats3$cell_count)
#head(df_gene_stats3)

write.csv(df_gene_stats3, paste0(results_dir, "/df_gene_stats_spleen_even.csv"), row.names = FALSE)



# BARPLOTS
head(df_gene_stats3)
table(df_gene_stats3$OfInterest_Reason)
# group by OfInterest_Reason to plot barplots across samples 
dim(df_gene_stats3)

# Step 1: Filter the data
filtered_df <- df_gene_stats3 %>%
  filter(gene %in% c("Cd3d", "Cd4", "Cd8a", "Cd68", "Cd163", "Cd19", "Itgax", "Itgam"))
  #filter(OfInterest_Reason == "Senescence") 

colnames(filtered_df)
dim(filtered_df)
filtered_df[,1]

# Step 2: Create the bar plot
ggplot(filtered_df, aes(x = gene, y = norm_expression[,1], fill = factor(sample))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Expression Percentage of Non-Zero Values for CellType Genes",
       x = "Gene",
       y = "Expression % Non-Zero",
       fill = "Sample") +
  theme_minimal()
ggsave(paste0(results_dir, "/barplot_celltype_genes.png"))

# for each unique OfInterest_Reason
for (reason in unique(df_gene_stats3$OfInterest_Reason)) {
  reason2 <- sanitize_name(reason)
  print(reason2)
  filtered_df <- df_gene_stats3 %>%
    filter(OfInterest_Reason == reason) 
  
  # Step 2: Create the bar plot
  ggplot(filtered_df, aes(x = gene, y = norm_expression[,1], fill = factor(sample))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = paste("Expression Percentage of Non-Zero Values for", reason2, "Genes"),
         x = "Gene",
         y = "Expression % Non-Zero",
         fill = "Sample") +
    theme_classic()
  ggsave(paste0(results_dir, "/barplot_", reason2, ".png"))
}


###################################
# plot grid of ImageDim Plots with select FOVs and gene combinations
###################################

# sobj_spleen1@images

# #sobj_spleen1@images$fov31 <- cropped.coords.fov
# Idents(sobj_spleen1) <- "Sample.ID"
# # The following is still not working; with subset fov, the molecules not plotted, or the entire slide is 
# plotted
# ImageDimPlot(sobj_spleen1, fov="fov_32", molecules = c("Mtor", "Cd3e")) 

ImageDimPlot(sobj_spleen5, fov="Sabaawynewcore091320245", axes=TRUE, flip_xy =FALSE )

#######################
# Define fovs, needed for spatial plots of molecules. 
# fov30 
# fov31, etc. 
df_meta <- sobj_spleen5@meta.data
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 30)
bbox_mtx
# BUG: The y and x coords are flipped in the Crop function
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], 
                       y = c(bbox_mtx[2,1], bbox_mtx[2,2]), x = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj_spleen5[["fov30"]] <- cropped.coords
# do we have the right subregion?
ImageDimPlot(sobj_spleen5, fov="fov30", axes=TRUE, flip_xy=FALSE)

# fov31
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 31)
bbox_mtx
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], 
                       y = c(bbox_mtx[2,1], bbox_mtx[2,2]), x = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj_spleen5[["fov31"]] <- cropped.coords
# fov32
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 32)
bbox_mtx
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], 
                       y = c(bbox_mtx[2,1], bbox_mtx[2,2]), x = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj_spleen5[["fov32"]] <- cropped.coords
# fov34
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 34)
bbox_mtx
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], 
                       y = c(bbox_mtx[2,1], bbox_mtx[2,2]), x = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj_spleen5[["fov34"]] <- cropped.coords
# fov35
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 35)
bbox_mtx
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], 
                       y = c(bbox_mtx[2,1], bbox_mtx[2,2]), x = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj_spleen5[["fov35"]] <- cropped.coords

DefaultBoundary(sobj_spleen5[["fov30"]]) <- "centroids"
DefaultBoundary(sobj_spleen5[["fov31"]]) <- "centroids" # "segmentation" # sets cell segmentation outline
DefaultBoundary(sobj_spleen5[["fov32"]]) <- "centroids" # sets cell segmentation outline
DefaultBoundary(sobj_spleen5[["fov34"]]) <- "centroids" # sets cell segmentation outline
DefaultBoundary(sobj_spleen5[["fov35"]]) <- "centroids"



#################
# Marker / Color lists
# 

genes_regions <- c("Cd3e", "Cd3d", "Cd3g", "Cd19", "Cd4", "Cd8a", "Cd8b1", "Cd68", "Cd163")
colors_regions <-  c("yellow", "yellow", "yellow","blue", "magenta", "magenta", "magenta", "orange", "orange")
names(colors_regions) <- genes_regions



genes_bcells 
colors_bcells_combined <- rep("blue", length(genes_bcells))
# Bcl2 is organge to highlight 
colors_bcells_sep <- c("blue", "lightblue", "green", "blue", "orange")
names(colors_bcells_combined) <- genes_bcells
names(colors_bcells_sep) <- genes_bcells
colors_bcells_combined
colors_bcells_sep


# endothelial
genes_vasc <- c("Pecam1", "Cd34")
colors_vasc <- c("red", "green")
names(colors_vasc) <- genes_vasc

# Create sample ggplot objects (you can replace these with your own plots)
# plot_list <- lapply(1:12, function(i) {
#   ImageDimPlot(sobj_spleen, fov="", label = TRUE)  
# })
Idents(sobj_spleen5) <- "Sample.ID"

plot_list <- list()

# Generate and store ggplots in the list

  ############################
  #  Region Markers 

  #molcol <-  c("yellow", "yellow", "yellow","blue")
  #names(molcol) <- c("Cd3e", "Cd3d", "Cd3g", "Cd19")
  p <- plot_molecules (sobj_spleen5, "fov30", mol_list = genes_regions, 
                       colors_regions, "Regions Sample 30 - WT", show_legend=FALSE)
  plot_list[[1]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_regions, 
                       colors_regions, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[2]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_regions, 
                       colors_regions, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[3]] <- ggplotGrob(p)
  

  ############################
  # B cells
  p <- plot_molecules (sobj_spleen5, "fov30", mol_list = genes_bcells, 
                       colors_bcells_combined, "B-cells  Sample 30 - WT", show_legend=FALSE)
  plot_list[[4]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_bcells, 
                       colors_bcells_combined, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[5]] <- ggplotGrob(p)
  
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_bcells, 
                       colors_bcells_combined, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[6]] <- ggplotGrob(p)
  
  ############################
  # Bcells separate colors
  p <- plot_molecules (sobj_spleen5, "fov30", mol_list = genes_bcells, 
                       colors_bcells_sep, "B cells (check Bcl2) Sample 30 - WT", show_legend=FALSE)
  plot_list[[7]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_bcells, 
                       colors_bcells_sep, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[8]] <- ggplotGrob(p)
  
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_bcells, 
                       colors_bcells_sep, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[9]] <- ggplotGrob(p)
  
  
  #############################
  # Vasc 
  
  p <- plot_molecules (sobj_spleen5, "fov30", mol_list = genes_vasc, 
                       colors_vasc, "Vasculature Sample 30 - WT", show_legend=FALSE)
  plot_list[[10]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_vasc, 
                       colors_vasc, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[11]] <- ggplotGrob(p)
  p
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_vasc, 
                       colors_vasc, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[12]] <- ggplotGrob(p)
  
  
length(plot_list)

# Convert ggplot objects to grobs
#grob_list <- lapply(plot_list, ggplotGrob)  

# Arrange the 12 plots in a 3x4 grid
# Open a PDF device to save the grid of plots
png(paste0(results_dir,"/grid_of_plots1.png"), width = 800, height = 1000)  # Specify the file name and dimensions
grid.arrange(grobs = plot_list, ncol = 3, nrow=4)
dev.off()

###################################
# Next set 
# Senescence page
genes_senescence[1:4] 
colors_senescence_sep <- c("blue", "red", "salmon", "magenta", "red", "orange", "blue")
names(colors_senescence_sep) <- genes_senescence


# IL/chemokines
genes_special_il
colors_il_combined <- rep("blue", length(genes_special_il))
colors_il_sep <- colorRampPalette(c("red", "blue", "green", "yellow","orange", "brown","grey20", "purple"))(length(genes_special_il))
names(colors_il_combined) <- genes_special_il
names(colors_il_sep) <- genes_special_il

# growth factors
genes_special_gf
colors_gf_combined <- rep("blue", length(genes_special_gf))
names(colors_gf_combined) <- genes_special_gf

plot_list2 <- list()

# Generate and store ggplots in the list

############################
#  Senescence SubTypes Markers 

#molcol <-  c("yellow", "yellow", "yellow","blue")
#names(molcol) <- c("Cd3e", "Cd3d", "Cd3g", "Cd19")
p <- plot_molecules (sobj_spleen5, "fov31", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Senescence I Sample 31 - WT", show_legend=FALSE)
plot_list2[[1]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[2]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[3]] <- ggplotGrob(p)

#####################
p <- plot_molecules (sobj_spleen5, "fov31", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Senescence II Sample 31 - WT", show_legend=FALSE)
plot_list2[[4]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[5]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[6]] <- ggplotGrob(p)


############################
# IL/chemokines
p <- plot_molecules (sobj_spleen5, "fov31", mol_list = genes_special_il, 
                     colors_il_combined, "SASP I Sample 31 - WT", show_legend=FALSE)
plot_list2[[7]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_special_il, 
                     colors_il_combined, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[8]] <- ggplotGrob(p)

# sobj, fov_name, mol_list, mol_colors, title
p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_special_il, 
                     colors_il_combined, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[9]] <- ggplotGrob(p)

############################
# IL/chemokines separate colors
p <- plot_molecules (sobj_spleen5, "fov31", mol_list = genes_special_il, 
                     colors_il_sep, "SASP I Sample 31 - WT", show_legend=FALSE)
plot_list2[[10]] <- ggplotGrob(p)

p <- plot_molecules (sobj_spleen5, "fov32", mol_list = genes_special_il, 
                     colors_il_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[11]] <- ggplotGrob(p)

# sobj, fov_name, mol_list, mol_colors, title
p <- plot_molecules (sobj_spleen5, "fov34", mol_list = genes_special_il, 
                     colors_il_sep, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[12]] <- ggplotGrob(p)


#############################
# growth factors 

# p <- plot_molecules (sobj, "fov31", mol_list = genes_special_gf, 
#                      colors_gf_combined, "Sample 31 - WT", show_legend=FALSE)
# plot_list[[10]] <- ggplotGrob(p)
# 
# p <- plot_molecules (sobj, "fov32", mol_list = genes_special_gf, 
#                      colors_gf_combined, "Sample 32 - Hetero", show_legend=FALSE)
# plot_list[[11]] <- ggplotGrob(p)
# 
# # sobj, fov_name, mol_list, mol_colors, title
# p <- plot_molecules (sobj, "fov34", mol_list = genes_special_gf, 
#                      colors_gf_combined, "Sample 34 -Homoz", show_legend=TRUE)
# plot_list[[12]] <- ggplotGrob(p)


length(plot_list2)

png(paste0(results_dir,"/grid_of_plots2.png"), width = 800, height = 1000) # , res=300 ) # , dpi=300)  # Specify the file name and dimensions
grid.arrange(grobs = plot_list2, ncol = 3, nrow=4)
dev.off()




############################
# Niches?
####### 
# Crop for slide 5, spleen only

colnames(sobj_spleen5@meta.data)
table(sobj_spleen5@meta.data$spatialclust_a57_1_assignments)

sobj_spleen5@images
table(sobj_spleen5@meta.data$spatialclust_af7_1_assignments)
Idents(sobj_spleen5) <- "spatialclust_af7_1_assignments"
#Idents(sobj_spleen) <- "spatialclust_a57_1_assignments" # not as good
# ImageDimPlot(sobj_spleen, fov = "spleenfov",cols = "polychrome",
#              coord.fixed = TRUE)

niche_plot_spleen(sobj_spleen5, "spleenfov",  "spatialclust_af7_1_assignments", "Niches in spleen")


## SpatialFeaturePlot
#sobj_test <- subset(sobj, subset = Sample.ID == 31))
cropped.coords <- Crop(sobj_spleen5[["Sabaawynewcore091320245"]], x = c(x_min, x_max), y = c(y_min, y_max), coords = "plot")

# add new fov
sobj_spleen5[["zoom"]] <- cropped.coords
sobj_spleen5@images$zoom 
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(sobj[["zoom"]]) <- "segmentation"
ImageDimPlot(sobj_spleen5, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Cd3e", "Cd3d"), nmols = 10000)




DefaultBoundary(sobj_spleen5[["fov30"]]) <- "segmentation" # sets cell segmentation outline
#DefaultBoundary(sobj[["fov30"]]) <- "tissue" # sets tissue outline
Idents(sobj) <- "Sample.ID"
molcol <-  c("yellow", "yellow", "yellow","blue")
names(molcol) <- c("Cd3e", "Cd3d", "Cd3g", "Cd19")
p <- ImageDimPlot(sobj_spleen5, fov = "fov31", axes = FALSE, # border.color = "white", border.size = 0.1, 
             cols = "polychrome",
             coord.fixed = TRUE, 
             molecules = c("Cd3e", "Cd3d", "Cd3g", "Cd19"), 
             nmols = 20000,
             mols.cols = molcol, 
             mols.size = 0.3, mols.alpha = 0.5,
             dark.background = FALSE) +
  ggtitle("Sample 31 - WT") +
  theme_minimal()
p


SpatialFeaturePlot(sobj_spleen5,  features = c("Cd3e", "Cd3d"), slot="counts") 
                   


##################################
# Boxplots
#
table(sobj_spleen@meta.data$Sample.Label)

plot_mouse_multiple_gene_boxplots <- function(df_gene_data, gene_list, group.by="Sample.ID"){
  
  if (group.by == "Sample.Label") {
    df_gene_data$Group <- factor(df_gene_data$Sample.Label, 
                                levels <- c("134 Spleen Normal",
                                 "158 Spleen Normal",
                                 "295 SpleenNormal",
                                 "296 Spleen Normal",
                                 "408 Spleen Normal",
                                 "680 Spleen Normal" 
                                 ))
  }
  else if (group.by == "Sample.ID") {
    df_gene_data$Group <- factor(df_gene_data$Sample.ID, 
                                 levels <- c("30",
                                             "31",
                                             "32",
                                             "33",
                                             "34",
                                             "35" 
                                 ))
  }
  else{
    df_gene_data$Group <- factor(df_gene_data[group.by])
  }
  
  Group <- group.by
  df_gene_data[[Group]] <- as.factor(df_gene_data[[Group]])
  # maker sure
  
  # filter out data not in groups
  df_gene_data <- df_gene_data %>% 
    dplyr::filter(Group %in% levels(Group))
  
  # Reshape the dataframe from wide to long format
  gene_data_long <- df_gene_data %>%
    pivot_longer(cols = gene_list, names_to = "Gene", values_to = "Expression")
  
  # Now you can plot using facet_wrap
  p <- ggplot(gene_data_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Gene, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() +
    labs(title = "Boxplots for Multiple Genes", x = "Group", y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = NA) 
    )
  return(p)
}

colnames(sobj_spleen@meta.data)

fetch_mouse_gene_data <- function(sobj_spleen, genes_list){
  # given a subset object and list of genes, fetch normalized data
  # and metadata needed for typical boxplots
  DefaultAssay(sobj_spleen) <- "RNA"
  norm_data <- FetchData(sobj_spleen, gene_names, layer = "counts") # s/b data
  df_meta <- sobj_spleen@meta.data[c('Sample.ID','Sample.Label', 'Patient')]
  head(df_meta)
  gene_data <- merge(norm_data, df_meta, how="inner", by="row.names")
  return (gene_data)
  
}
Layers(sobj_spleen)
DefaultAssay(sobj_spleen) <- "SCT"

# Squamous # too low expression to tell 
file_suffix <- "Tcells"
gene_names <- c("Cd3d","Cd3e", "Cd3g", "Mki67" ) # TP63
gene_data1 <- fetch_mouse_gene_data(sobj_spleen, gene_names)
head(gene_data1)
p <- plot_mouse_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle(paste0("T Cell markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Tcells.png"), p, dpi = 300)

# repeat for combos
file_suffix <- "Senescence"
gene_names <- genes_senescence
gene_data1 <- fetch_mouse_gene_data(sobj_spleen, gene_names)
head(gene_data1)
p <- plot_mouse_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle("Senescense Cell markers in spleen")
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,".png"), p, dpi = 300)



# histogram of cell counts for gene Mki67
# norm:
hist(gene_data1$Mki67)
length(gene_data1[gene_data1$Mki67 > 0, ])
length(gene_data1$Mki67)

gene_counts <- FetchData(sobj_spleen, gene_names, layer="counts")
head(gene_counts)
dim(gene_counts[gene_counts$Mki67 > 0, ])

gene_norm <- FetchData(sobj_spleen, gene_names, layer="data")
dim(gene_norm[gene_norm$Mki67 > 0, ])
gene_norm[gene_norm$Mki67 > 0, ]

# we *DO* have 13K positive values. 
table(gene_norm$Mki67)


# Calculate the mean and the quartiles
mean_val <- mean(gene_norm$Mki67)
std_dev <- sd(gene_norm$Mki67)
#quartiles <- quantile(gene_norm$Mki67, probs = c(0.25, 0.75))

mean_val
std_dev

# Plot the box plot with whiskers showing the quartiles
# Create some sample data
set.seed(123)
data <- data.frame(value = rnorm(100))


# Calculate mean and standard deviation
mean_val <- mean(data$value)
std_dev <- sd(data$value)
std_dev

# Plot showing mean and mean ± 1 standard deviation
ggplot(data, aes(x = "", y = value)) +
  geom_boxplot(aes(ymin = mean_val - std_dev, lower = mean_val - std_dev, 
                   middle = mean_val, upper = mean_val + std_dev, 
                   ymax = mean_val + std_dev), stat = "identity") +
  annotate("point", x = 1, y = mean_val, color = "red", size = 3) +  # Highlight mean
  labs(title = "Box Plot with Mean and ±1 Standard Deviation", 
       x = "", y = "Value") +
  theme_minimal()

######## AGain
MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

MeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

g1 <- ggplot(mtcars, aes(factor(am), mpg)) + geom_boxplot() +
  ggtitle("Regular Boxplot")

g2 <- ggplot(gene_data1, aes(factor(Sample.ID), Mki67)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="red") + 
  ggtitle("Boxplot: Min, Mean-1SEM, Mean, Mean+1SEM, Max")

MinMeanSEMMax(gene_data1$Mki67)

g3 <- ggplot(gene_data1, aes(factor(Sample.ID), Mki67)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="red") + 
  ggtitle("Boxplot: Mean-1SEM, Mean, Mean+1SEM")


# Plot showing mean and mean ± 1 standard deviation
ggplot(gene_norm, aes(x = "", y = Mki67)) +
  geom_boxplot(aes(ymin = mean_val - std_dev, lower = mean_val - std_dev, 
                   middle = mean_val, upper = mean_val + std_dev, 
                   ymax = mean_val + std_dev), stat = "identity") +
  geom_point(aes(x = "", y = mean_val), color = "red", size = 3) +  # Highlight mean
  labs(title = "Box Plot with Mean and ±1 Standard Deviation", 
       x = "", y = "Value") +
  theme_minimal()

# Create a density plot showing only the non-zero values
ggplot(gene_counts, aes(x = Mki67)) +
  geom_density() +
  xlim(0, NA) +  # Focus on non-zero values
  labs(title = "Density Plot of Non-Zero Values", x = "Value", y = "Density")

#################################################
# InSituType Cell Typing
# See the vignette code here for more details:
# https://github.com/Nanostring-Biostats/InSituType
#################################################
#library(InSituType)
#data("ioprofiles")
#data("iocolors")
#data("mini_nsclc")
set.seed(0)



############################
# InSituType Cell Typing
# caution: loading the biomaRt package will break InSituType by overriding 'select'
# depends on normalized Seurat object, and vector of negative probe mean expresison per cell
# 
############################

#readRDS(paste0(results_dir, "/sobj_from_Load_v1.RDS"))
sobj123 <- readRDS(paste0(saved_objs, "/sobj123_w_InSituType_UMAP.RDS")) 

numeric_vec <- as.numeric(as.character(negmeans_vec2$count))
names(numeric_vec) <- rownames(negmeans_vec2)
negmeans_vec <- numeric_vec
negmeans_vec[1:5]
sum(negmeans_vec) # 4002
#insitu_matrix_Az <- as.matrix(read.csv( "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_AzimuthLung6_12.csv",
#                                        row.names = 1))


#############################
# optional additional barplots by celltype, e.g. for different grouping 
df_meta <- sobj@meta.data
df_meta$cell_type <- df_meta$NS_Insitu_Celltype.level1
group.by <- "Run_Tissue_name"

# group by tissue and then cell_type 
df_meta <- df_meta %>%
  dplyr::mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
  dplyr::select(!!group.by, cell_type)
#         tissue  cell_type
#c_1_1_1    CPA    tumor.9
#c_1_1_2    CPA neutrophil
#c_1_1_3    CPA   tumor.13

df_celltype <- df_meta %>%
  group_by_at(c(group.by, "cell_type")) %>%
  summarise(count = n(), 
            #log_count = log2(n() + 1), 
            .groups = 'drop')

table(df_celltype$cell_type)
my_colors = c("Tumor" = "red", "Immune" = "blue", "Fibroblast" = "green", "Endothelial" = "purple", "Epithelial" = "orange")
p <- plot_celltypes(df_celltype, group.by, "InSituType", colors=my_colors, log2=FALSE)
p
ggsave(paste0(results_dir, "/celltype_barplot_by_", group.by, ".png"), p, width = 10, height = 5, dpi = 300)
##############################

# Cell counts by Sample
group.by <- "PatientID"
df_meta <- sobj@meta.data
df_meta$cell_type <- df_meta$NS_Insitu_Celltype.level1
# group by tissue and then cell_type 
df_meta <- df_meta %>%
  dplyr::mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
  dplyr::select(!!group.by, cell_type)


df_celltype <- df_meta %>%
  group_by_at(c(group.by, "cell_type")) %>%
  summarise(count = n(), .groups = 'drop')

plot_celltypes(df_celltype, group.by, "ALL")


df_meta <- sobj@meta.data
plot_grouped_barplot(df_meta, cell.type='NS_Insitu_Celltype.level1', group.by='PatientID', title=" ")


# Save the Seurat object
saveRDS(sobj, file = paste0(saved_objs, "/sobj_w_InSituType.RDS"))
#sobj <- readRDS(paste0(results_dir, "/sobj_w_InSituType.RDS"))


insitu_matrix_cosMx_IO <- as.matrix(read.csv( "/Users/annmstrange/Documents/git/CosMx-Cell-Profiles/Human/IO/IO.profiles.csv",
                                           row.names = 1))
head(insitu_matrix_cosMx_IO)

# 
# insitu_matrix_safeTME <- safeTME 
dim(group_averages_matrix)
colnames(group_averages_matrix) <- "Tumor"
head(group_averages_matrix)
# dim(safeTME)
# # how many genes are in common?
common_genes <- intersect(rownames(group_averages_matrix), rownames(insitu_matrix_cosMx_IO))
length(common_genes) # 839

# join mtx and safeTME by common genes, keeping all mtx genes
#mtx_common <- mtx[common_genes,]
IO_common <- insitu_matrix_cosMx_IO[common_genes,]


########################################
# PCA, UMap etc 
# 1. altogether vs 
# 2. patient cells only

########################################

DefaultAssay(sobj) <- "RNA"
# 
VariableFeatures(sobj) <- c(VariableFeatures(obj_list[[1]]), 
                            VariableFeatures(obj_list[[2]]),
                            VariableFeatures(obj_list[[3]]), 
                            VariableFeatures(obj_list[[4]]))

# assume FindVariableFeatures has been run
p1 <- VariableFeaturePlot(sobj)
top10 <- head(VariableFeatures(sobj))
p2 <- LabelPoints(plot = p1, points=top10, repel=TRUE)
p1 + p2
## Try each slide separately
p1 <- VariableFeaturePlot(obj_list[[1]])
top10 <- head(VariableFeatures(obj_list[[1]]))
p2 <- LabelPoints(plot = p1, points=top10, repel=TRUE)
p1 + p2



sobj <- RunPCA(sobj, verbose = FALSE, npcs = 20)

# Assess how PCA did with different npcs, 20, 30, 50 e.g. 
# ElbowPlot to decide best cutoff.  DimPlot colored by possible batch effect variables. e.g. slides should look integrated.
ElbowPlot(sobj)


