###########################
## ---------------------------
##
## Script name: mouse_spleen_RNA_EDA.R
##
## Purpose of script: explore the Seurat data object produced by the AtoMx pipeline,
##   generate plots, evaluate rna quality controls
##.  summarize and display what has already been done in AtoMx
## 1. Use Seurat LoadNanostring() method (vs exported Seurat RDS) because it includes FOVs 
##    1b. subtract Negative Probes from the counts matrix
##    1c. Apply QC filters to the Seurat object
## 2. Load Sample metadata per cell:
##    2a. FOV and Sample.ID assignments
##    2b. Apply panCK/CD45 phenotypes to cells where possible. 
##    2c. overlay the pathology annotations to get additional meta: path_annot_id and name
## 3. Apply Cell Typing to the Seurat object with InSituType 
## 4. Split sobj into tissue, patient, and samples vs Organoids
## 5. Run SCTransform, PCA, UMAP on the Seurat object(s)
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

old.wd <- getwd()           # save the current working directory
setwd("~/Documents/git/cosMx-tools/")   # set working directory (mac)
source("R/seurat_utils.R") # load up the functions we need
source("R/image_reg_tools.R")
#source("R/gene_annotations.R") # needs a little fixing
source("R/cell_typing_utils.R")

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
rna_dir <- paste0(root_root_dir, '/', slide_name, '/', run_name, '/AnalysisResults/s2wi20fpqr') 
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
df_gene_data <- get_nanostring_gene_annotations(ns_gene_fn)

# drop ns_annotations
#df_gene_data <- df_gene_data %>% dplyr::select(-ns_annotations)

df_special_genes <- read.csv(paste0(markers_dir, '/GenesOfInterest_RNA_HES.csv'))

colnames(df_special_genes) <- c('Display_Name', 'OfInterest_Reason', "Notes", "GroupAs", "Highlight")
df_special_genes <- df_special_genes %>% dplyr::select(Display_Name, OfInterest_Reason, Notes, GroupAs, Highlight)
head(df_special_genes)

genes_special <- df_special_genes$Display_Name

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
sobj_load <- readRDS(paste0(saved_rds_dir, "/latest_sobj_15Oct24.rds"))
sobj <- sobj_load

##################################
# Marker lists
#################################
#iocolors

heatmap_markers = "GLF3"
lc_markers = c("EGFR", "TP53")
immune_markers <- c("PTPRC", "CD3E", "CD3D", "CD3G", "CD4", "CD8A", "CD8B", "FOXP3")
epithelial_markers <- c("KRT5", "KRT6A", "KRT6B", "KRT13", "TP63")
#adeno_markers <- c("MUC1", "MUC5AC", "MUC5B", "KRT7", "KRT20", 'ROS1', 'CLDN3', 'NKX2-1')
#adeno_markers_gentles <- c("MUC1", "MUC5B", "KRT7", "KRT20", 'ROS1', 'CLDN3', 'NKX2-1')
basal_keratins <- c("KRT5", "KRT6A", "KRT6B", "KRT13", "KRT14", "TP63")

# lineage markers from out list (matched to RNA 1K)
lineage_endothelial <- c("PECAM1", "CD34")
lineage_fibroblast <- c("FN1")
lineage_macrophage <- c("CD68", "CD163")
lineage_macrophage_polarization <- c("CD8A", "PDL1", "PD1",  "INOS")
lineage_epithelial <- c("EPCAM", "KRT5", "KRT6A", "KRT6B", "KRT13", "TP63")

gentles_need_alias <- c('NTS',	'GPX2',	'CSTA',	'DAPL1',	'VSNL1','DSG3',	'SERPINB1',	'DSC3',	'CLCA2',	'MIR205HG',	'TP63',	'PKP1',			'TRIM29',	'AKR1C2',	'PTHLH',		'NT5E',		'ACSL5',	'NDNF',		'ARSE',	'CCNJL',	'GOLT1A',	'MUC1',	'CAPN8',	'MLPH',			'SLC44A4',	'LMO3',	'FOLR1','NAPSA','SFTA2',	'NKX2-1',	'SFTA3',	'TOX3',	'CLDN3','TESC',	'CRIP2')
# TODO: use biomaRt to get alternate possible names for these
# KRT6A/B/C has a problem w / 
gentles_adeno <- c('GPNMB',	'KRT13',	'KRT14',	'KRT16',	'KRT5',	'NTRK2',	'S100A8',	'SERPINB5')
gentles_sq <- c('ADGRF5',	'AGR2',	'AREG',	'CEACAM6',	'DPP4')





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

celltype_colors <- c("B.cell"= "darkblue", 
              
              "Fibroblast" = "#999999" ,
              "Macrophage" =  "#006600" ,
              "Mast" = "#FFFF00" ,
              "mDC" = "#00FF00"  ,
              "Monocyte" = "#33CC00"  ,
              "Neutrophil" =  '#B3DE69', 
              "NK" =    '#80B1D3', 
              "pDC" = "#00FFFF" ,
              "Plasmablast" =   "#3399CC"  ,
              "T.cell" = 'skyblue',
              "T.CD4.memory" =   'purple',
              "T.CD4.naive" =   '#BC80BD',
              "T.CD8.memory" =     "#996633" ,
              "T.CD8.naive" =      "pink",
              "Treg" =    "#FF66FF",
              "Endothelial"= '#FDB462',   
              "Tumor" = "#FF0000",
              "Stroma" ='#8DD3C7',
              "a" = "grey90",
              "b" = 'wheat',
              'c' = 'peachpuff1',
              'd' = '#FF9900',
              'e' = 'yellow',
              'f' = 'grey40',
              'g' = 'grey30',
              'h' = 'grey20',
              'i' = 'peachpuff3',
              'j' = 'pink',
              "k" = "grey90",
              "l" = 'grey80',
              'm' = 'grey70',
              'n' = 'grey60',
              'o' = 'grey50',
              'p' = '#FF0000',
              'q' = 'grey30',
              'r' = 'grey20',
              's' = 'grey10',
              't' = 'pink'
              )


known_color_values = c("B_lineage"='#8DD3C7', 
                       "DC"='yellow', 
                       "eo/baso/mast"='purple', "erythroid"='#FB8072', "HSC/MPP and pro"='#80B1D3', 
                       "MK"='#FDB462', "monocyte"='#B3DE69', "neutrophil"='magenta', 
                       "stroma"='blue', "T/NK"='#BC80BD')


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
  
  obj_list[[i]] <- LoadNanostring(data.dir = flatfiles_dir, fov = slide_names[i], assay="Nanostring")
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
  obj_list[[i]] = remove_sys_control(obj_list[[i]])
  print(paste0("Num Features after removing System Controls (expect 960): ", nrow(obj_list[[i]])))
  
  #obj_list[[i]][["Nanostring"]]@counts <- remove_sys_control(obj_list[[i]], "Nanostring", "counts")
  
  # Negative controls should be subtracted/removed before normalization, but not before running InSituType
  # However, we keep 0 values for Negative controls in the counts matrix after subtracting, and set aside the means for later
  negmeans_list[[i]] <- get_neg_control_means(obj_list[[i]])
  # add metadata for neg
  obj_list[[i]] <- AddMetaData(obj_list[[i]], metadata = data.frame(neg = negmeans_list[[i]]))
 
  #SetAssayData(obj_list[[i]], layer = "counts", new.data = subtract_neg_control(obj_list[[i]], "Nanostring", "counts"))
  obj_list[[i]] <- subtract_neg_control(obj_list[[i]])
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
  obj_list[[i]]
  
}

# Merge Seurat objects using chaining
# Merging allows for group metadata to be added to the Seurat object, and also allows for
# joint dimensional redux and clustering on the underlying RNA expression data if desired. (but we won't)
sobj <- merge(obj_list[[1]], y = c(obj_list[[2]]))
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
#sobj <- readRDS(paste0(results_dir, "/sobj_from_Load_v1.RDS"))


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
df_fov_meta <- df_sample[, c('Run_Tissue_name', 'fov', 'Sample.ID')]
  

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

calculate_mouse_PanCK_PT <- function(value) {
  # if (is.na(value) || value > 4000) {
  #   return('na')
  if (value > 1200) {
    return(1)
  } else {
    return(0)
  }
}

calculate_mouse_CD45_PT <- function(value) {
  if (is.na(value)) {
    return(0)
  } else
    if (value > 3000) {
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
add_mouse_sample_metadata <- function(obj, df_fov_meta){
  
  df_metadata <- obj@meta.data
  # Apply the function to create the PanCK.PT column
  
  
  df_metadata2 <- merge(df_metadata, df_fov_meta, 
                        by.x = c('fov', 'Run_Tissue_name'), by.y = c('fov', 'Run_Tissue_name'),
                        suffixes = c(".x",""), how = 'inner')
  print(head(df_metadata2))
  
  # df_metadata2$Mean.PanCK[is.na(df_metadata2$Mean.PanCK)] <- 0
  # df_metadata2$Max.PanCK[is.na(df_metadata2$Max.PanCK)] <- 0
  # df_metadata <- df_metadata %>%
  #  mutate(PanCK.PT = sapply(Mean.PanCK, calculate_PanCK_PT))
  
  # sample thresholds
  panCK.thresholds <- c("30" = 1200, "31" = 1200, "32" = 1200, "33" = 1200, "34" = 1200, "34" = 1200,
                        "35" = 1200, "0" = 1200)
  df_metadata2 <- df_metadata2 %>%
    dplyr::mutate(
      PanCK.threshold = panCK.thresholds[Sample.ID],  # Look up threshold by Sample.ID
      PanCK.PT = as.integer(Mean.PanCK >= PanCK.threshold)  # Threshold comparison
    )
  print(table(df_metadata2$PanCK.PT))
  df_metadata2$Mean.CD45[is.na(df_metadata2$Mean.CD45)] <- 0
  df_metadata2$Max.CD45[is.na(df_metadata2$Max.CD45)] <- 0
  CD45.thresholds <-  c("30" = 3000, "31" = 3000, "32" = 3000, "33" = 3000, "34" = 3000, "34" = 3000,
  "35" = 3000, "0" = 3000)
  
  df_metadata <- df_metadata %>%
   mutate(CD45.PT = sapply(Mean.CD45, calculate_CD45_PT))
  df_metadata2 <- df_metadata2 %>%
    mutate(CD45.threshold = CD45.thresholds[Sample.ID],          # Look up threshold by sample
           CD45.PT = as.integer(Mean.CD45 >= CD45.threshold))
  
  print(table(df_metadata2$CD45.PT))

  cols_to_add <- c("cell_id", "CD45.PT", "PanCK.PT", "Mean.PanCK", "Mean.CD45",
                  "PanCK.threshold", "CD45.threshold", colnames(df_fov_meta))
  
  cols_to_add <- c("cell_id", colnames(df_fov_meta))
  
  print(paste("cells before merge", nrow(df_metadata2)))
  
  #df_metadata2$PanCK.PT[is.na(df_metadata2$PanCK.PT)] <- 1
  #df_metadata2$CD45.PT[is.na(df_metadata2$CD45.PT)] <- 0
  
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
  print(length(cell_IDs_not_in_metadata)) # expect 0 
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


sobj <- add_mouse_sample_metadata(sobj, df_fov_meta2)
table(sobj@meta.data$Organ)
table(sobj@meta.data$Genotype)

# Exclude FOVs marked for exclusion in the fov metadata
table(sobj@meta.data$Exclude) # 2823 cells to remove
sobj <- subset(sobj, subset = Exclude != "x")


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
ImageDimPlot(sobj, fov="Sabaawynewcore091320246")

Idents(sobj) <- "CD45.PT"
ImageDimPlot(sobj, fov="Sabaawynewcore091320246")
ImageDimPlot(sobj, fov="Sabaawynewcore091320245")

Idents(sobj_spleen) <- "Genotype"
# Todo: use subset of cells for performance
ImageDimPlot(sobj_spleen, fov="Sabaawynewcore091320246", col=c("red","blue","green"))

table(sobj_spleen@meta.data$Sample.ID)
sobj_spleen1 <- subset(sobj_spleen, subset = Sample.ID %in% c(31,32,34) 
                       & Run_Tissue_name == "Sabaawy new core 09/13/2024 5")

Idents(sobj_spleen1) <- "Genotype"
ImageDimPlot(sobj_spleen1, fov="Sabaawynewcore091320245")
# FeaturePlot
ImageFeaturePlot(sobj_spleen1, fov="Sabaawynewcore091320245",
                 features = c("Vim", "Cd3e"),
                 dark.background = FALSE)




#############################
# Sample Level Stats
#############################

# TODO: move this below path annotations
df_sample_stats <- sobj@meta.data %>%
  group_by(Sample.ID) %>%  # Group by Sample.ID
  summarise(
    cell_count = n(),               # Count the number of rows per Sample.ID
    unique_fovs = n_distinct(fov),   # Count the number of unique fovs
    PanCK_tumor_count = sum(PanCK.PT == 1),
    Path_annot_tumor_count = sum(path_annot_name.1 == "Annotated Tumor"),
    #PanCK_pct_tumor = (PanCK_tumor_count / cell_count) * 100  # Percentage of tumor cells (PanCK.PT == 1)
  )
head(df_sample_stats)


sobj_spleen <- subset(sobj, subset = Organ == "Spleen" & Sample.ID >= "30")
table(sobj_spleen@meta.data$Sample.ID)
# Do we have any SCTransforms?
sobj_spleen@assays$SCT@SCTModel.list 

# normalize for spleen only
sobj_spleen <- SCTransform(sobj_spleen, assay = "Nanostring", verbose = TRUE)

saveRDS(sobj_spleen, paste0(saved_rds_dir,"/sobj_spleen.rds"))
# get counts of each molecule by sample.

exp_data <- as.matrix(GetAssayData(sobj_spleen, assay="SCT", layer="counts"))
norm_exp_data <- as.matrix(GetAssayData(sobj_spleen, assay="SCT", layer="data"))
class(exp_data)
dim(exp_data) # 1000 x 191K
# how to subset matrix by gene

# which genes in genes_special not found as rownames in exp_data
genes_special_not_found <- setdiff(genes_special, rownames(exp_data))
print(genes_special_not_found)

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

df_meta <- sobj_spleen@meta.data
# get cell counts per sample
cell_counts_by_sample <- table(df_meta$Sample.ID)

pct_non_zero <- function(x) {
  length(x[x >0]) / length(x) * 100
}

mean_non_zero <- function(x) {
  mean(x[x >0])
}


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

write.csv(df_gene_stats3, paste0(results_dir, "/df_gene_stats_allspleen.csv"), row.names = FALSE)

###################################
# plot grid of ImageDim Plots
###################################
library(gridExtra)

##########################
# first, make some sub-fovs

# fov_nm e.g. "Sabaawynewcore091320245"
get_coords_fov <- function(sobj, sample_ids, fov_nm) {
  df_meta <- sobj@meta.data
  df_meta2 <- df_meta %>%
    filter(Sample.ID %in% sample_ids)
  
  padding <- 0
  x_min <- round(min(df_meta2$CenterX_global_px)) - padding
  y_min <- round(min(df_meta2$CenterY_global_px)) - padding
  x_max <- round(max(df_meta2$CenterX_global_px)) + padding
  y_max <- round(max(df_meta2$CenterY_global_px)) + padding
  print(paste(y_min, y_max, x_min, x_max))
  
  sobj_fov <- subset(sobj, subset = Sample.ID %in% sample_ids )
  
  cropped.coords.fov <- Crop(sobj[[fov_nm]], 
                             x = c(x_min, x_max), y = c(y_min,y_max), coords = "tissue")
  
  # Returns an fov object which should be added to the original object. 
  return(cropped.coords.fov)
}

newfov <- get_coords_fov(sobj_spleen, c(30,31,32,34,35), "Sabaawynewcore091320245")
sobj_spleen@images$spleenfov <- newfov

# fov <- get_coords_fov(sobj_spleen1, 31, "Sabaawynewcore091320245")
# sobj_spleen1$fov_31 <- fov
# 
# fov <- get_coords_fov(sobj_spleen1, 32, "Sabaawynewcore091320245")
# sobj_spleen1$fov_32 <- fov
# 
# fov <- get_coords_fov(sobj_spleen1, 34, "Sabaawynewcore091320245")
# sobj_spleen1$fov_34 <- fov
# 
# 
# sobj_spleen1@images

# 
# #sobj_spleen1@images$fov31 <- cropped.coords.fov
# Idents(sobj_spleen1) <- "Sample.ID"
# # still not working; with subset fov, the molecules not plotted.
# ImageDimPlot(sobj_spleen1, fov="fov_32", molecules = c("Mtor", "Cd3e")) 


#######################
# Define fovs
# fov31
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 31)
bbox_mtx
cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], 
                       x = c(bbox_mtx[2,1], bbox_mtx[2,2]), y = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj[["fov31"]] <- cropped.coords
# fov32
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 32)
bbox_mtx
cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], 
                       x = c(bbox_mtx[2,1], bbox_mtx[2,2]), y = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj[["fov32"]] <- cropped.coords
# fov34
bbox_mtx <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 34)
bbox_mtx
cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], 
                       x = c(bbox_mtx[2,1], bbox_mtx[2,2]), y = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj[["fov34"]] <- cropped.coords

DefaultBoundary(sobj[["fov31"]]) <- "centroids" # "segmentation" # sets cell segmentation outline
DefaultBoundary(sobj[["fov32"]]) <- "centroids" # sets cell segmentation outline
DefaultBoundary(sobj[["fov34"]]) <- "centroids" # sets cell segmentation outline


sobj6 <- subset (sobj, subset = Sample.ID %in% c(6))
df_meta6 <- sobj6@meta.data
bbox_mtx <- get_bbox_of_sample(df_meta6, "Sabaawy new core 09/13/2024 5", 6)
bbox_mtx
cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], 
                       x = c(bbox_mtx[2,1], bbox_mtx[2,2]), y = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
sobj[["fov6"]] <- cropped.coords
DefaultBoundary(sobj[["fov6"]]) <- "centroids"
genes_regions6 <- c("Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd8b1", "Foxp3")
colors_regions6 <-  c("red", "red", "red", "green", "green", "blue")
names(colors_regions6) <- genes_regions6

p <- plot_molecules (sobj, "fov6", mol_list = genes_regions6, 
                colors_regions6, "Regions Sample 6", show_legend=TRUE, molsize=1)
p
ggsave(paste0(results_dir,"regions_sample6.png"), p, width=10, height=10)

# cell counts by fov
count_mtx6 <- FetchData(sobj, fov = "fov6", vars=genes_regions6, layer = "counts")
head(count_mtx6)

file_suffix <- "Brain6"
gene_data1 <- fetch_mouse_gene_data(sobj6, gene_names6)
head(gene_data1)
p <- plot_mouse_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle("Sample 6")
p

DefaultAssay(sobj6) <- "Nanostring"
counts_data <- FetchData(sobj6, genes_regions6, layer = "counts") # s/b data
df_meta6 <- sobj6@meta.data[c('Sample.ID','Sample.Label', 'Patient', 'fov')]
head(df_meta6)
gene_data <- merge(counts_data, df_meta6, how="inner", by="row.names")
plot_mouse_multiple_gene_boxplots(gene_data, genes_regions6)+
  ggtitle("Sample 6")

head(gene_data)

# I want to plot the number of rows with Cd3e > 0, Cd3d > 0, Cd3g > 0, or Cd8a > 0 
# for each fov
# "Cd8a", "Cd8b1"

summary_df <- gene_data %>%
  filter(Foxp3 > 0 ) %>%  # Filter rows where any of the conditions are true
  group_by(fov) %>%  # Group by fov
  summarise(FoxP3_pos_cells = n())  

summary_df



#################
# Marker / Color lists
# 

genes_regions <- c("Cd3e", "Cd3d", "Cd3g", "Cd19", "Cd4", "Cd8a", "Cd8b1", "Cd68", "Cd163")
colors_regions <-  c("yellow", "yellow", "yellow","blue", "magenta", "magenta", "magenta", "orange", "orange")
names(colors_regions) <- genes_regions



genes_bcells 
colors_bcells_combined <- rep("blue", length(genes_bcells))
# Bcl2 is organge to highlight 
colors_bcells_sep <- c("blue", "lightblue", "green", "orange")
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
Idents(sobj) <- "Sample.ID"

plot_list <- list()

# Generate and store ggplots in the list

  ############################
  #  Region Markers 

  #molcol <-  c("yellow", "yellow", "yellow","blue")
  #names(molcol) <- c("Cd3e", "Cd3d", "Cd3g", "Cd19")
  p <- plot_molecules (sobj, "fov31", mol_list = genes_regions, 
                       colors_regions, "Regions Sample 31 - WT", show_legend=FALSE)
  plot_list[[1]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj, "fov32", mol_list = genes_regions, 
                       colors_regions, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[2]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj, "fov34", mol_list = genes_regions, 
                       colors_regions, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[3]] <- ggplotGrob(p)
  

  ############################
  # B cells
  p <- plot_molecules (sobj, "fov31", mol_list = genes_bcells, 
                       colors_bcells_combined, "B-cells  Sample 31 - WT", show_legend=FALSE)
  plot_list[[4]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj, "fov32", mol_list = genes_bcells, 
                       colors_bcells_combined, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[5]] <- ggplotGrob(p)
  
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj, "fov34", mol_list = genes_bcells, 
                       colors_bcells_combined, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[6]] <- ggplotGrob(p)
  
  ############################
  # Bcells separate colors
  p <- plot_molecules (sobj, "fov31", mol_list = genes_bcells, 
                       colors_bcells_sep, "B cells (check Bcl2) Sample 31 - WT", show_legend=FALSE)
  plot_list[[7]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj, "fov32", mol_list = genes_bcells, 
                       colors_bcells_sep, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[8]] <- ggplotGrob(p)
  
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj, "fov34", mol_list = genes_bcells, 
                       colors_bcells_sep, "Sample 34 -Homoz", show_legend=TRUE)
  plot_list[[9]] <- ggplotGrob(p)
  
  
  #############################
  # Vasc 
  
  p <- plot_molecules (sobj, "fov31", mol_list = genes_vasc, 
                       colors_vasc, "Vasculature Sample 31 - WT", show_legend=FALSE)
  plot_list[[10]] <- ggplotGrob(p)
  
  p <- plot_molecules (sobj, "fov32", mol_list = genes_vasc, 
                       colors_vasc, "Sample 32 - Hetero", show_legend=FALSE)
  plot_list[[11]] <- ggplotGrob(p)
  
  # sobj, fov_name, mol_list, mol_colors, title
  p <- plot_molecules (sobj, "fov34", mol_list = genes_vasc, 
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
p <- plot_molecules (sobj, "fov31", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Senescence I Sample 31 - WT", show_legend=FALSE)
plot_list2[[1]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov32", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[2]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov34", mol_list = genes_senescence[1:4] , 
                     colors_senescence_sep, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[3]] <- ggplotGrob(p)

#####################
p <- plot_molecules (sobj, "fov31", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Senescence II Sample 31 - WT", show_legend=FALSE)
plot_list2[[4]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov32", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[5]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov34", mol_list = genes_senescence[5:7] , 
                     colors_senescence_sep, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[6]] <- ggplotGrob(p)


############################
# IL/chemokines
p <- plot_molecules (sobj, "fov31", mol_list = genes_special_il, 
                     colors_il_combined, "SASP I Sample 31 - WT", show_legend=FALSE)
plot_list2[[7]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov32", mol_list = genes_special_il, 
                     colors_il_combined, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[8]] <- ggplotGrob(p)

# sobj, fov_name, mol_list, mol_colors, title
p <- plot_molecules (sobj, "fov34", mol_list = genes_special_il, 
                     colors_il_combined, "Sample 34 -Homoz", show_legend=TRUE)
plot_list2[[9]] <- ggplotGrob(p)

############################
# IL/chemokines separate colors
p <- plot_molecules (sobj, "fov31", mol_list = genes_special_il, 
                     colors_il_sep, "SASP I Sample 31 - WT", show_legend=FALSE)
plot_list2[[10]] <- ggplotGrob(p)

p <- plot_molecules (sobj, "fov32", mol_list = genes_special_il, 
                     colors_il_sep, "Sample 32 - Hetero", show_legend=FALSE)
plot_list2[[11]] <- ggplotGrob(p)

# sobj, fov_name, mol_list, mol_colors, title
p <- plot_molecules (sobj, "fov34", mol_list = genes_special_il, 
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

colnames(sobj_spleen@meta.data)
table(sobj_spleen@meta.data$spatialclust_a57_1_assignments)

sobj_spleen@images
table(sobj_spleen@meta.data$spatialclust_af7_1_assignments)
Idents(sobj_spleen) <- "spatialclust_af7_1_assignments"
#Idents(sobj_spleen) <- "spatialclust_a57_1_assignments" # not as good
# ImageDimPlot(sobj_spleen, fov = "spleenfov",cols = "polychrome",
#              coord.fixed = TRUE)

niche_plot_spleen(sobj_spleen, "spleenfov",  "spatialclust_af7_1_assignments", "Niches in spleen")


## SpatialFeaturePlot
#sobj_test <- subset(sobj, subset = Sample.ID == 31))
cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], x = c(x_min, x_max), y = c(y_min, y_max), coords = "plot")

# add new fov
sobj[["zoom"]] <- cropped.coords
sobj@images$zoom 
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(sobj[["zoom"]]) <- "segmentation"
ImageDimPlot(sobj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Cd3e", "Cd3d"), nmols = 10000)




DefaultBoundary(sobj[["fov30"]]) <- "segmentation" # sets cell segmentation outline
#DefaultBoundary(sobj[["fov30"]]) <- "tissue" # sets tissue outline
Idents(sobj) <- "Sample.ID"
molcol <-  c("yellow", "yellow", "yellow","blue")
names(molcol) <- c("Cd3e", "Cd3d", "Cd3g", "Cd19")
p <- ImageDimPlot(sobj, fov = "fov31", axes = FALSE, # border.color = "white", border.size = 0.1, 
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


SpatialFeaturePlot(sobj,  features = c("Cd3e", "Cd3d"), slot="counts") 
                   


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
  DefaultAssay(sobj_spleen) <- "SCT"
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
#negmeans_vec2 <- readRDS(paste0(results_dir, "/negmeans_vec2.RDS")) 
#readRDS( paste0(results_dir, "/merged_mtx.RDS")) # merged_mtx

# bug in here.. 
# negmean <- unlist(negmeans_list[1:4])
# negmean[1:5]
# negmeans_vec[1:5]
# # IF negmens_list is not available as a variable use the following to load from file. 
# negmeans_vec2 <- read.table(paste0(results_dir, "/negmeans_list.txt"), row.names=1, col.names = "count")
# # ignore warning
#class(negmeans_vec2)
numeric_vec <- as.numeric(as.character(negmeans_vec2$count))
names(numeric_vec) <- rownames(negmeans_vec2)
negmeans_vec <- numeric_vec
negmeans_vec[1:5]
sum(negmeans_vec) # 4002
#insitu_matrix_Az <- as.matrix(read.csv( "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_AzimuthLung6_12.csv",
#                                        row.names = 1))

# from here Supplementary Tables: S10 the InsituType matrix used by the authors for NSCLC samples
# https://www.nature.com/articles/s41587-022-01483-z#Sec30
#insitu_matrix_He <- as.matrix(read.csv( "/Volumes/T7Shield/cosMxGeoMx_Analysis/reference_data/NSCLC_CellProfileMatrix.csv",
#                                        row.names=1))

#head(insitu_matrix_He, row.names=1)
#insitu_matrix <- insitu_matrix_He
#rownames(insitu_matrix)
#colnames(insitu_matrix)
#insitu_matrix_He <- insitu_matrix_He[, colnames(insitu_matrix_Az)]
#head(insitu_matrix)

#rownames(insitu_matrix_Az[1:10,])
#dim(insitu_matrix_Az)

# Trying different matrices
# safeTME focuses on immune cells, HCA healthy lung cells, CPA lung cancer cell lines 
# What He 2022 clearly did was identify cells that were certainly each tumor cell type, and then used the
# methods to create an InSituType matrix to classify the cells in the NSCLC samples. (like ID'ing anchor cells)
# 
#insitu_matrix_safeTME <- as.matrix(read.csv( "/Volumes/T7Touch4/public_data/CellProfileLibrary/Human/Adult/ImmuneTumor_safeTME_profilematrix.csv",
#                                             row.names = 1))
# insitu_matrix_HCA <- as.matrix(read.csv( "/Volumes/T7Touch4/public_data/CellProfileLibrary/Human/Adult/Lung_HCA_profilematrix.csv",
#                                          row.names = 1))
# 
# insitu_matrix_CPA <- as.matrix(read.csv( "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_CPA_LC_Lines.csv",
#                                           row.names = 1))

# join the matrices
#rownames(insitu_matrix_safeTME)[1:10]
#rownames(insitu_matrix_CPA)[1:10]
# outer join by rownames

# sup <- run_insitu_type(sobj, negmean,title_suffix="All", # colname="cell_type", 
#                 insitu_matrix=insitu_matrix, results_dir=results_dir)

# All
# possible error: "unable to find an inherited method for function ‘select’ for signature ‘"data.frame"’"
# try loading this package again?: 

# InSituType has Dependency on SCTransform for normalization and on negmean for the negative mean values
#sobj <- SCTransform(sobj,assay = "Nanostring", verbose = FALSE)
# Supervised Cell Typing 

# O.G. with He 2022 matrix
# sum(negmeans_vec) 
#insitu_matrix <- insitu_matrix_He

cohort_psr01 <- subset(sobj123, subset = fov == "PSR01")

cohort_data1 <- sobj@meta.data[, colnames(sobj@meta.data) %>%
                                grep ("Mean.PanCK|Mean.CD45|PanCK.PT|CD45.PT|Area.um2", ., value = TRUE)] %>%
                                filter(sobj@meta.data$Run_Tissue_name == 'PSR-01')

cohort_data2 <- sobj@meta.data[, colnames(sobj@meta.data) %>%
                                 grep ("Mean.PanCK|Mean.CD45|PanCK.PT|CD45.PT|Area.um2", ., value = TRUE)] %>%
  filter(sobj@meta.data$Run_Tissue_name == 'PSR-02')

cohort_data3 <- sobj@meta.data[, colnames(sobj@meta.data) %>%
                                 grep ("Mean.PanCK|Mean.CD45|PanCK.PT|CD45.PT|Area.um2", ., value = TRUE)] %>%
  filter(!sobj@meta.data$Run_Tissue_name %in% c('PSR-10','PSR-02'))
#cohort_data <- sobj@meta.data[, c("PanCK.PT", "CD45.PT")]
cohort_data3[1:5, ]
cohort_data3[,'Mean.PanCK'] <- 1000
cohort_data3[,'PanCK.PT'] <- 1
cohort_data3[,'Mean.CD45'] <- 80


# convert NAs to defaults for tumor cells
cohort_data[is.na(cohort_data$Mean.PanCK), "PanCK.PT"] <- 1000
cohort_data[is.na(cohort_data$Mean.CD45), "CD45.PT"] <- 0
cohort_data[is.na(cohort_data$PanCK.PT), "PanCK.PT"] <- 1
cohort_data[is.na(cohort_data$CD45.PT), "CD45.PT"] <- 0
cohort_data[is.na(cohort_data$Area.um2), "Area.um2"] <- 0
#cohort_data$PT <- cohort_data$PanCK.PT - cohort_data$CD45.PT
table(cohort_data$PanCK.PT)
# perform automatic cohorting
#cohort <- fastCohorting(cohort_data[,c("Mean.PanCK", "Mean.CD45")],
#                        gaussian_transform = TRUE)

cohort1 <- fastCohorting(as.matrix(cohort_data1[,c("Mean.PanCK", "Mean.CD45", "Area.um2")]), 
                         gaussian_transform = TRUE)
cohort2 <- fastCohorting(as.matrix(cohort_data2[,c("Mean.PanCK", "Mean.CD45", "Area.um2")]), 
                         gaussian_transform = TRUE)
cohort3 <- fastCohorting(as.matrix(cohort_data2[,c("Mean.PanCK", "Mean.CD45", "Area.um2")]), 
                         gaussian_transform = TRUE)
cohort_data <- rbind(cohort_data1, cohort_data2, cohort_data3)
cohort_data[1:5,]

table(cohort2)
cohort <- fastCohorting(as.matrix(cohort_data1[,c("Mean.PanCK", "Mean.CD45", "Area.um2")]), 
                         gaussian_transform = TRUE)
# Error when more than one slide are involved, likely due to too many identical values? 
#Error in mclust::predict.Mclust(object = mc, newdata = mat) : 
#  object not of class 'Mclust'

# ("Gaussian_transform = TRUE" maps variables to gaussians in order to 
#  place dramatically different variables on the same scale.)
table(cohort) # 3 cohorts shows possibly good separation

# 
# 
# 
# # Generalized cell types:
# df_meta_to_add <- df_meta_to_add %>%
#   mutate(
#     NS_Insitu_Celltype.level2 = case_when(
#       NS_Insitu_Celltype %in% c("T.cell.CD4", "T.cell.CD8", "T.cell.regulatory","Treg", "T.CD8.naive", "T.CD8.memory", "T.CD4.memory", "T.CD4.naive") ~ "T.cell",
#       NS_Insitu_Celltype %in% c("B.cell", "plasmablast", "Plasma", "Plasmacytoid.dendritic.cell", "Dendritic.cell","Plasmablast") ~ "B.cell",
#       NS_Insitu_Celltype %in% c("macrophage", "Macrophage", "monocyte", "Monocyte") ~ "Macrophage",
#       NS_Insitu_Celltype %in% c("NK", "Neutrophil","neutrophil", "NK.cell","pDC", "mDC", "Mast.cell","mast") ~ "Innate immune",
#       NS_Insitu_Celltype %in% c("epithelial") ~ "Epithelial",
#       NS_Insitu_Celltype %in% c("fibroblast") ~ "Fibroblast",
#       NS_Insitu_Celltype %in% c("endothelial") ~ "Endothelial",
#       NS_Insitu_Celltype %in% c("tumor.6") ~ "SqCC",
#       NS_Insitu_Celltype %in% c("Tumor") ~ "Tumor",
#       NS_Insitu_Celltype %in% c("tumor.9", "tumor.12", "tumor.13", "tumor.5") ~ "Adeno",
#       TRUE ~ NS_Insitu_Celltype
#     )
#   )
# Generalized cell types:
# df_meta_to_add <- df_meta_to_add %>%
#   mutate(
#     NS_Insitu_Celltype.level1 = case_when(
#       NS_Insitu_Celltype.level2 %in% c("T.cell","B.cell", "macrophage", "innate immune") ~ "immune",
#       #NS_Insitu_Celltype.level2 %in% c("epithelial") ~ "epithelial",
#       NS_Insitu_Celltype.level2 %in% c("fibroblast") ~ "fibroblast",
#       NS_Insitu_Celltype.level2 %in% c("tumor", "adeno", "SqCC", "epithelial") ~ "tumor",
#       TRUE ~ NS_Insitu_Celltype
#     )
#   )
# df_meta_to_add <- df_meta_to_add %>%
#   mutate(
#     NS_Insitu_Celltype.level1 = case_when(
#       NS_Insitu_Celltype.level2 %in% c("T.cell","B.cell", "Macrophage", "Innate immune") ~ "Immune",
#       #NS_Insitu_Celltype.level2 %in% c("epithelial") ~ "epithelial",
#       NS_Insitu_Celltype.level2 %in% c("Fibroblast") ~ "Fibroblast",
#       NS_Insitu_Celltype.level2 %in% c("Tumor", "Adeno", "SqCC", "epithelial") ~ "Tumor",
#       TRUE ~ NS_Insitu_Celltype
#     )
#   )
# 
# table(df_meta_to_add$NS_Insitu_Celltype.level1)
# 
# sobj <- AddMetaData(sobj, df_meta_to_add)


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

plot_grouped_barplot <- function(df_meta, cell.type='NS_Insitu_Celltype.level1', group.by='Run_Tissue_name', title=""){
  df_meta$cell_type <- df_meta[[cell.type]]
  
  df_meta <- df_meta %>%
    dplyr::mutate(!!group.by := gsub(" ", "_", .data[[group.by]])) %>%
    dplyr::select(!!group.by, cell_type)
  
  df_celltype <- df_meta %>%
    group_by_at(c(group.by, "cell_type")) %>%
    summarise(count = n(), .groups = 'drop')
  
  p <- plot_celltypes(df_celltype, group.by, title, log2=FALSE)
  #ggsave(paste0(results_dir, "/celltype_barplot_by_", group.by, ".png"), p, width = 10, height = 5, dpi = 300)
  
  return(p)
}
df_meta <- sobj@meta.data
plot_grouped_barplot(df_meta, cell.type='NS_Insitu_Celltype.level1', group.by='PatientID', title=" ")


# Save the Seurat object
saveRDS(sobj, file = paste0(saved_objs, "/sobj_w_InSituType.RDS"))
#sobj <- readRDS(paste0(results_dir, "/sobj_w_InSituType.RDS"))

#######################################
# Get only anchor tumor cells InSituType
######################################

sobj12 <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02"))
# narrow down cells to only tumor cells
sobj1 <- subset(sobj12, subset = fov %in% c(24,27,22,17,4))

# liver fov 17 has much higher EPCAM, use higher TH. 

feature_names <- rownames(sobj1[["SCT"]])
counts_matrix <- sobj1[["SCT"]]@counts
counts_matrix[1:10, 1:5]

# cells where EpCAM >= 1
counts_matrix[feature_names == "EPCAM", 1:5]

epcam_cells <- which(counts_matrix[feature_names == "EPCAM", ] >= 1)
length(epcam_cells) # 2588

# which are also PTPRC and CD31 negative
ptprc_cells <- which(counts_matrix[feature_names == "PTPRC", ] >= 1)
length(ptprc_cells) # 617

epcam_cells2 <- which(counts_matrix[feature_names == "EPCAM", ] >= 1 &
                        counts_matrix[feature_names == "PTPRC", ] == 0 
)

epcam_cells2_liver <- which(counts_matrix[feature_names == "EPCAM", ] >= 2 &
                        counts_matrix[feature_names == "PTPRC", ] == 0 
)
                        #counts_matrix[feature_names == "CD31", ] == 0)
length(epcam_cells2) # 2376

sobj_tumor <- subset(sobj1, subset = (cell_id %in% names(epcam_cells2) & fov %in% c(24,27,22,4)) |
                      (cell_id %in% names(epcam_cells2_liver) & fov == 17 ))

dim(sobj_tumor)
sobj_tumor@meta.data$cell_id[1:5]
table(sobj_tumor@meta.data$Sample.Label)


Idents(sobj_tumor) <- "path_annot_name"
# p_immune_molecules <- ImageDimPlot(sobj_tumor, 
#                                    fov="PSR01",
#                                    cells = row.names(s_lung@meta.data)[which(s_lung@meta.data$fov == 24)],
#                                    alpha = 0.5,
#                                    molecules = c("PTPRC", "CD8A", "CD4", "CD68", "CD19", "FOXP3", "Mean.CD45"), #
#                                    #"EPCAM", "KRT5", "KRT7", "Mean.PanCK")
#                                    mols.size = 0.02,
#                                    axes = FALSE) + 
#   ggtitle("Immune Molecules in Lung - Bug: showing full slide, not just fov")
# p_immune_molecules

# ImageDimPlot
sobj_tumor_22 <- subset(sobj_tumor, subset = fov == 22)
ImageDimPlot(sobj_tumor_22, 
             fov = "PSR01",  
             size = 4,
            )

fovs_to_filter = c(32,69,90,161,117,141,202,201,193,6)
fovs_for_cellid = c(168,210,179,180)
sobj2 <- subset(sobj12, subset = fov %in% c(fovs_to_filter,fovs_for_cellid) & Run_Tissue_name == "PSR-02")

sobj2@meta.data$cell_id[1:5] # "c_4_6_7" cell Id format
# for specific cell ids.
# liver (4_1)
# fov 168: get cells 17,25,27,31,43,46,309,345
# 
# hemorhaged area:  5_2 pleura_rx
# fov 210, get cells 190,192,193
# fov 179, get cells 49,66, 92,147,276
# fov 180, get cells 1,3,8,49,96,121

feature_names <- rownames(sobj2[["SCT"]])
counts_matrix <- sobj2[["SCT"]]@counts

epcam_cells <- which(counts_matrix[feature_names == "EPCAM", ] >= 1)
length(epcam_cells) # 2588

# which are also PTPRC and CD31 negative
ptprc_cells <- which(counts_matrix[feature_names == "PTPRC", ] >= 1)
length(ptprc_cells) # 617

epcam_cells2 <- which(counts_matrix[feature_names == "EPCAM", ] >= 1 &
                        counts_matrix[feature_names == "PTPRC", ] == 0 
)
#counts_matrix[feature_names == "CD31", ] == 0)
epcam_names = names(epcam_cells2)
length(epcam_cells2) # 332

cells168 <- paste0("c_4_168_", c(17,25,27,31,43,46,309,345))
cells210 <- paste0("c_4_210_", c(190,192,193))
cells179 <- paste0("c_4_179_", c(49,66, 92,147,276))
cells180 <- paste0("c_4_180_", c(1,3,8,49,96,121))
tumor_anchor_cells <- c(cells168,cells210,cells179,cells180, epcam_names)
length(tumor_anchor_cells) # 1506
tumor_anchor_cells # 
sobj2_tumor <- subset(sobj2, subset = cell_id %in% tumor_anchor_cells) # 364

df_meta2 <- sobj2_tumor@meta.data
table(df_meta2$Sample.Label)

sobj2_tumor_179 <- subset(sobj2_tumor, subset = fov == 179)
ImageDimPlot(sobj2_tumor_179, 
             fov = "PSR02",  
             size = 4,
)



cell_type_list <- sobj_tumor$orig.ident #sobj_tumor$Sample.ID
mtx <- sobj_tumor@assays$SCT@counts
mtx2 <- sobj2_tumor@assays$SCT@counts
head(mtx)
mtx3 <- cbind(mtx, mtx2)
mtx <- mtx3

group_averages <- calculate_group_averages(mtx, cell_type_list)
group_averages_matrix <- as.matrix(group_averages)

#cell_type_list <- sobj2_tumor$orig.ident   #sobj2_tumor$Sample.ID
#mtx2 <- sobj2_tumor@assays$SCT@counts
#group_averages2 <- calculate_group_averages(mtx, cell_type_list) # detailed

# tumor group only
#group_averages2 <- calculate_group_averages(mtx, cell_type_list)
#group_averages_matrix2 <- as.matrix(group_averages2)

# this only works if theres no combined columns 
# mtx <- cbind(group_averages_matrix, group_averages_matrix2)

# then cosMx IO
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

#######################
# Merge Cell Profile matrices
#######################

merged_mtx <- merge_matrices(group_averages_matrix, insitu_matrix_cosMx_IO)


#mtx_safeTME <- cbind(mtx, safeTME_common)
#dim(mtx_safeTME)
# Convert row names to a column
# mtx1_df <- data.frame(Gene = rownames(group_averages_matrix), group_averages_matrix, row.names = NULL)
# mtx2_df <- data.frame(Gene = rownames(IO_common), IO_common, row.names = NULL)
# # Perform the left outer join
# merged_df <- merge(mtx1_df, mtx2_df, by = "Gene", all.x = TRUE)
# # Convert back to matrix if needed (excluding the 'Gene' column)
# merged_mtx <- as.matrix(merged_df[, -1])
# # Set the row names to be the 'Gene' column
# rownames(merged_mtx) <- merged_df$Gene
# 
# # fill any NA
# merged_mtx[is.na(merged_mtx)] <- 0

dim(merged_mtx)
# Now, merged_mtx contains all rows from mtx1 matched with mtx2

write.csv(as.data.frame(merged_mtx), 
          "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_Sample_Tumor1.csv",
          row.names = TRUE)


df_metadata <- sobj@meta.data
table(df_metadata$Run_Tissue_name)
sobj123 <- subset(sobj, subset = Run_Tissue_name != "PSR-CPA")
counts <- GetAssayData(sobj123, slot = "data", layer="SCT") %>%
  as.matrix() %>%
  t()

# negmean needs same rows as counts
rownames(counts)[1:5]
negmean <- negmean[rownames(counts)]

dim(counts)
length(negmean)
# error: caused by a rowsum being NA.  avoid this. 
any(rowSums(counts) == 0) # NA values in counts are a problem 
# unsup <- insitutypeML(x = counts,
#                     neg = negmean,
#                     cohort = NULL, # cohort_data,
#                     reference_profiles = merged_mtx)   
# 86 genes in the count data are missing from reference_profiles and will be omitted from cell typing (ok-ish)
# Error: Not compatible with requested type: [type=list; target=double]



########################################
# PCA, UMap etc 
# 1. altogether vs 
# 2. patient cells only

########################################

DefaultAssay(sobj) <- "SCT"
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
sobj123 <- subset(sobj, subset = Run_Tissue_name %in% c('PSR-01', 'PSR-02', 'PSR-Organoids'))
sobj123 <- RunPCA(sobj123, verbose = FALSE, npcs = 20)
ElbowPlot(sobj123)

# s_lung <- subset(sobj, subset = Run_Tissue_name %in% c('PSR-01') & tissue == 'lung')
# s_lung_24 <- subset(s_lung, subset = Run_Tissue_name %in% c('PSR-01') & tissue == 'lung' & fov %in% c('24'))
# s_lung_27 <- subset(s_lung, subset = Run_Tissue_name %in% c('PSR-01') & tissue == 'lung' & fov %in% c('27'))
# 
# 
# ElbowPlot(s_lung)
# ImageDimPlot(s_lung, fov="PSR01",  cols = "glasbey", dark.background = FALSE)
# ImageDimPlot(s_lung_24, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE)
# ggsave(paste0(results_dir, "ImageDimPlot_lung_fov24_4Clusters.png"), last_plot())
# ImageDimPlot(s_lung_27, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE)
# ggsave(paste0(results_dir, "ImageDimPlot_lung_fov27_4Clusters.png"), last_plot())
# 
# # same plots with markers
# Idents(s_lung_24) <- 'EPCAM'
# ImageDimPlot(s_lung_24, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE) # , features = c('EPCAM', 'PTPRC', 'CD31', 'CD10'))
# 
# ImageFeaturePlot(s_lung_24, features="EPCAM", size=2) # fov warning but worked 
# ImageFeaturePlot(s_lung_24, features="EPCAM", size=2, dark.background = FALSE) # fov warning but worked 
# ImageFeaturePlot(s_lung_24, fov="PSR01", features="EPCAM", size=2, dark.background = FALSE)
# ggsave(paste0(results_dir, "ImageFeaturePlot_lung_fov24_EPCAM.png"), last_plot())

# something went wrong with plotting too much here
# Warning: No FOV associated with assay 'SCT', using global default FOV
#ImageDimPlot(s_lung_24, molecules="EPCAM", size=5, mols.cols = "red", dark.background = FALSE)

DimPlot(sobj, reduction = "umap")



# min_dist 0.1, n_neighbors doesn't make a differene


# lists (imprecise)
lymphocyte_cells <- c('CD8A', 'CD8B', 'CD4', 'CD3E', 'CD3G', 'CD3D')
ctyotoxic <- c('GZMA', 'GZMB', 'PRF1', 'NKG7', 'CCL5','CTSW', 'GNLY','KLRB1','KLRD1','KLRK1','CCL13','CD290')
macrophage <- c('CD68', 'CD14', 'CD163', 'CD86', 'CD206', 'CD11b', 'CD11c', 'CD64', 'CD16', 'CD32', 'CD163', 'CD169', 'CD204', 'CD206')
mast <- c('MS4A2', 'TPSAB1', 'CPA3', 'HDC', 'TPSB2')
neutrophil <- c('CD66b', 'CD11b', 'CD16', 'CD32', 'CD64', 'CD66b', 'CD11b', 'CD16', 'CD32', 'CD64')
nk <- c('NCR1', 'NCR2', 'NCR3', 'NKp30', 'NKp44', 'NKp46', 'CD56', 'CD57', 'CD94', 'CD158', 'CD159', 'CD314', 'CD336', 'CD337', 'CD355', 'CD336', 'CD337', 'CD355')
bcell <- c('CD19', 'CD20', 'CD22', 'CD24', 'CD27', 'CD38', 'CD40', 'CD45', 'CD79A', 'CD79B', 'CD80', 'CD81', 'CD83', 'CD86', 'CD138', 'CD180', 'CD185', 'CD196', 'CD197', 'CD200', 'CD220', 'CD221', 'CD222', 'CD223', 'CD224', 'CD225', 'CD226', 'CD227', 'CD228', 'CD229', 'CD230', 'CD231', 'CD232', 'CD233', 'CD234', 'CD235', 'CD236', 'CD237', 'CD238', 'CD239', 'CD240', 'CD241', 'CD242', 'CD243', 'CD244', 'CD245', 'CD246', 'CD247', 'CD248', 'CD249', 'CD250', 'CD251', 'CD252', 'CD253', 'CD254', 'CD255', 'CD256', 'CD257', 'CD258', 'CD259', 'CD260', 'CD261', 'CD262', 'CD263', 'CD264', 'CD265', 'CD266', 'CD267', 'CD268', 'CD269', 'CD270', 'CD271', 'CD272', 'CD273', 'CD274', 'CD275', 'CD276', 'CD277', 'CD278', 'CD279', 'CD280', 'CD281', 'CD282', 'CD283', 'CD284', 'CD285', 'CD286', 'CD287', 'CD288', 'CD289', 'CD290', 'CD291', 'CD292', 'CD293', 'CD294', 'CD295', 'CD296', 'CD297', 'CD298', 'CD299', 'CD300', 'CD301', 'CD302', 'CD303', 'CD304', 'CD305', 'CD306', 'CD307', 'CD308', 'CD309', 'CD310', 'CD311', 'CD312', 'CD313', 'CD314', 'CD315', 'CD316', 'CD317', 'CD318', 'CD319', 'CD320', 'CD321', 'CD322')
endothelial <- c('CD31', 'CD34', 'CD105', 'CD144', 'CD146')
epithelial <- c('CD326', 'CD324', 'CD325', 'CD326', 'CD327', 'CD328', 'CD329', 'CD330', 'CD331', 'CD332', 'CD333', 'CD334', 'CD335', 'CD336', 'CD337', 'CD338', 'CD339', 'CD340', 'CD341', 'CD342', 'CD343', 'CD344', 'CD345', 'CD346', 'CD347', 'CD348', 'CD349', 'CD350', 'CD351', 'CD352', 'CD353', 'CD354', 'CD355', 'CD356', 'CD357', 'CD358', 'CD359', 'CD360', 'CD361', 'CD362', 'CD363', 'CD364', 'CD365', 'CD366', 'CD367', 'CD368', 'CD369', 'CD370', 'CD371', 'CD372', 'CD373', 'CD374', 'CD375', 'CD376', 'CD377', 'CD378', 'CD379', 'CD380', 'CD381', 'CD382', 'CD383', 'CD384', 'CD385', 'CD386', 'CD387', 'CD388', 'CD389', 'CD390', 'CD391', 'CD392', 'CD393', 'CD394', 'CD395', 'CD396', 'CD397', 'CD398', 'CD399', 'CD400', 'CD401', 'CD402', 'CD403', 'CD404', 'CD405', 'CD406', 'CD407', 'CD408', 'CD409', 'CD410', 'CD411', 'CD412', 'CD413', 'CD414', 'CD415', 'CD416', 'CD417', 'CD418', 'CD419', 'CD420', 'CD421', 'CD422', 'CD423', 'CD424', 'CD425', 'CD426', 'CD427', 'CD428', 'CD429', 'CD430', 'CD431', 'CD432', 'CD433', 'CD434', 'CD435', 'CD436', 'CD437', 'CD438', 'CD439', 'CD440', 'CD441', 'CD442', 'CD443', 'CD444', 'CD445', 'CD446', 'CD447')
fibroblasts <- c('CD90', 'CD105', 'CD140a', 'CD140b', 'CD146', 'CD166', 'COL1A1', 'COL3A1', 'COL6A1', 'COL6A2', 'DCN','GREM1', 'PAMR1', 'TAGLN')
tcell <- c('CD3', 'CD4', 'CD8', 'CD25', 'CD27', 'CD28', 'CD45', 'CD45RA', 'CD45RO', 'CD69', 'CD127', 'CD152', 'CD154', 'CD195', 'CD196', 'CD197', 'CD198', 'CD199', 'CD200', 'CD201', 'CD202', 'CD203', 'CD204', 'CD205', 'CD206', 'CD207', 'CD208', 'CD209', 'CD210', 'CD211', 'CD212', 'CD213', 'CD214', 'CD215', 'CD216', 'CD217', 'CD218', 'CD219', 'CD220', 'CD221', 'CD222', 'CD223', 'CD224', 'CD225', 'CD226', 'CD227', 'CD228', 'CD229', 'CD230', 'CD231', 'CD232', 'CD233', 'CD234', 'CD235', 'CD236', 'CD237', 'CD238', 'CD239', 'CD240', 'CD241', 'CD242', 'CD243', 'CD244', 'CD245', 'CD246', 'CD247', 'CD248', 'CD249', 'CD250', 'CD251', 'CD252', 'CD253', 'CD254', 'CD255', 'CD256', 'CD257', 'CD258', 'CD259', 'CD260', 'CD261', 'CD262', 'CD263', 'CD264', 'CD265', 'CD266', 'CD267', 'CD268', 'CD269', 'CD270', 'CD271', 'CD272', 'CD273', 'CD274', 'CD275', 'CD276', 'CD277', 'CD278', 'CD279', 'CD280', 'CD281', 'CD282', 'CD283', 'CD284', 'CD285', 'CD286', 'CD287', 'CD288', 'CD289', 'CD290', 'CD291', 'CD292', 'CD293', 'CD294', 'CD295', 'CD296', 'CD297', 'CD298', 'CD299', 'CD300', 'CD301', 'CD302', 'CD303', 'CD304')


adeno <- c('KRT7','KRT20','TP53','RB1')
mor_genes <- c('MET', 'EGFR', 'TP53','ERBB2','ERBB3')
scc_genes <- c('TP63', 'KRT5', 'KRT14', 'KRT17', 'KRT6A/B/C')
sclc_genes <- c('NCAM1','FGFR1','MYC')
  
# Question: in paired samples, did these go up or down?




# group markers into aggregate features
mtx <- GetAssayData(s_lung, layer = "data")[lymphocyte_cells,,drop=FALSE ]
class(mtx)

# subscript out of bounds means a marker wasn't found
# Error: No cell overlap between new meta data and Seurat object means we're summing on wrong axis.
#s_lung[["lymph"]] <- colSums(GetAssayData(s_lung, layer="data")[lymphocyte_cells,,drop=FALSE])






DoHeatmap(sobj, features=TopFeatures(sobj, dim=1, nfeatures=20))
clusters <- Idents(s_lung)
clusters
VlnPlot(s_lung, features = c("PTPRC", "CD3E"))
VlnPlot(s_lung, features = c("PTPRC", "CD3E",))

FeaturePlot(s_lung, features = c("PTPRC", "CD3E", "CD4", "CD8A", "CD19"))
# "CD14", "CD56", "CD16", "CD11c", "CD11b", "CD123", "CD141", "CD1c)
FeaturePlot(s_lung, features = c("MGP", "DCN", "LUM", "EPCAM", "CD68", "SLPI"))


# subsample to make plot readable 
sampled_cells <- sample(colnames(sobj123), 10000)


subset123 <- subset(sobj123, cells = sampled_cells)
DimPlot(subset123, reduction = "pca", group.by = "Run_Tissue_name", alpha=0.5, cols='glasbey')
ggsave(paste0(results_dir, "/PCA_plot_Run_Tissue_Name123_20pcas.png"), plot = last_plot(), width = 10, height = 10, dpi = 300)


# Look for the 3 sliaes bascially overlapping; do not want separation by slide
#Idents(subset123) <- "Sample.ID"
DimPlot(subset123, reduction = "pca", group.by = "Sample.ID", alpha=0.5, cols='glasbey')
ggsave(paste0(results_dir, "/PCA_plot_SampleID123_20pcas.png"), plot = last_plot(), width = 10, height = 10, dpi = 300)
# same plot without CPA

# see with Cell Pellet Arrays (more distinctions; ok)
#DimPlot(sobj, reduction = "pca", group.by = "Sample.ID")
DimPlot(subset123, reduction = "pca", group.by = "tissue", alpha=0.5, cols='glasbey')

# Lung only, tumor cells (by panCK only)
subset_lung <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") & tissue == "lung" & Mean.PanCK >= 100)
dim(subset_lung)

sobj@images
ImageDimPlot(subset_lung, fov = "PSR01", group.by = "Sample.ID")
ImageDimPlot(subset_lung, group.by = "Sample.ID", cols = c("Mean.PanCK"), ncol = 3, nrow = 1)

#subset_Organoids <- subset(sobj, subset = Run_Tissue_name == "PSR-Organoids")
DimPlot(subset_lung, reduction = "pca", group.by = "Sample.ID", alpha=0.5, cols='polychrome')

sobj <- FindNeighbors(sobj, dims = 1:20)
sobj <- FindClusters(sobj,  cluster.name='seurat_cluster', verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:20)

p <- DimPlot(sobj, reduction="umap", group.by= c("ident", "Run_Tissue_name"), alpha=0.5, cols='polychrome')
show(p)
ggsave(paste0(results_dir, "Umap_plot1.png"), plot = p, width = 10, height = 10, dpi = 300)

p <- DimPlot(sobj, reduction="umap", group.by= c("spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"))
show(p)
ggsave(paste0(results_dir, "Umap_plot_by_niche.png"), plot = p, width = 10, height = 10, dpi = 300)

# UMAP plot with only patient cells
sobj123 <- subset(sobj, subset = Run_Tissue_name %in% c('PSR-01', 'PSR-02', 'PSR-Organoids'))
DimPlot(sobj123, reduction="umap", group.by= c("ident", "Run_Tissue_name"), alpha=0.5, cols='polychrome')
DimPlot(sobj123, reduction="umap", group.by= c("tissue", "Sample.ID"), alpha=0.5, cols='polychrome')
DimPlot(sobj123, reduction="umap", group.by= c("tissue", "Sample.ID"), alpha=0.5, cols='polychrome')
# Why is T2 vertebrae adding so much heterogeneity? 

# Again, tumor only
sobj_tumor <- subset(sobj, subset = Run_Tissue_name %in% c('PSR-01', 'PSR-02') & Mean.PanCK >= 100)
DimPlot(sobj_tumor, reduction="umap", group.by= c("tissue", "Sample.ID"))

# Lung only, tumor cells (by panCK only)
subset_lung <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") & tissue == "lung" & Mean.PanCK >= 50)
dim(subset_lung)
DimPlot(sobj_lung, reduction="umap", group.by= c("ident", "Run_Tissue_name"))


cell_names <- Cells(sobj)[1:1000]
cell_names
DimPlot(sobj, reduction="umap", cells = cell_names, split.by=c("Run_Tissue_name"))

Idents(sobj) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
DimPlot(sobj, reduction="umap", split.by=c("Run_Tissue_name"))



##############################
# InSituType w cosMx IO semisupervised
###############################3

# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/vignette-basic-analysis/assets/4.-cell-typing.html
refprofiles <- read.csv("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv", row.names = 1, header = TRUE)
refprofiles <- refprofiles[is.element(rownames(refprofiles), colnames(counts)), ]
p <- pheatmap::pheatmap(sweep(refprofiles, 1, pmax(apply(refprofiles, 1, max), 0.2), "/"), 
                   col = colorRampPalette(c("white", "darkblue"))(100))

p 
ggsave(paste0(results_dir, "/cosMxIO_refprofiles_heatmap.png"), plot = p, width = 10, height = 10, dpi = 300)


# do I want this plot?
for (cell in unique(res$clust)) {
  png(paste0(mydir, "/results/celltypev1_", cell, ".png"), 
      width = diff(range(xy[,1]))*.7, height = diff(range(xy[,2]))*.7, units = "in", 
      res = 400)  # res of 400 is pretty good; 600 is publication-quality
  par(mar = c(0,0,0,0))
  plot(xy, pch = 16, col = scales::alpha(cols[res$clust], 0.3), cex = 0.1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(xy[res$clust == cell, ], pch = 16, cex = 0.1, col = "black")
  legend("top", legend = cell)
  dev.off()
}


# Cell typing deserves careful attention: it’s a complex process, and your entire analysis depends on its results. Scrutinize the above QC plots for suspicious patterns. These could be:
#   
#   Unlikely spatial distributions
# Unlikely cell type abundances, e.g. a far too many immune cells, or not enough of a basic cell type
# Implausible marker genes expressed at high levels
# Biologically dissimilar cell types with similar expression profiles (near each other in UMAP space, confused with each other in the flightpath, showing similar profiles in the heatmap.)
# For a deep dive on this topic, see the Insitutype FAQs.

# 
# 
# ##############################
# # InSituType w cosMx IO PLUS our tumors semisupervised
# ###############################
# 
# #saveRDS(sobj, file = paste0(results_dir, "/sobj_w_InSituType2.RDS"))
# 
sobj123 <- subset(sobj, subset = Run_Tissue_name != "PSR-CPA")
# we already ran SCTransform on each slide
#sobj123 <- SCTransform(sobj123, assay = "Nanostring", verbose = TRUE)

sobj123 <- PrepSCTFindMarkers(sobj123, assay = "SCT", verbose = TRUE)
DefaultAssay(sobj123) <- "SCT"
sobj123 <- FindVariableFeatures(sobj123, selection.method = "vst", nfeatures = 1000) # all
# Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
p1a <- VariableFeaturePlot(sobj123) + ggtitle("All patient cells")
top10 <- head(VariableFeatures(s_lung))
p2a <- LabelPoints(plot = p1a, points=top10, repel=TRUE)
p2a
ggsave(paste0(results_dir, "VariableFeaturePlot_patientTissue.png"), last_plot())

# DefaultAssay(sobj123) <- "Nanostring"
#sobj123 <- SCTransform(sobj123, assay = "Nanostring", verbose = TRUE)
sobj123 <- RunPCA(sobj123, verbose = FALSE, npcs = 50)
ElbowPlot(sobj123) # 15 or 20 PCS is a good cutoff
sobj123 <- RunPCA(sobj123, verbose = FALSE, npcs = 20)

# run this to correct for 3 separate SCTransforms
sobj123 <- PrepSCTFindMarkers(sobj123, assay = "SCT", verbose = TRUE)
#s_tumor_2 <- SCTransform(s_tumor_2, assay="Nanostring", verbose = FALSE)

sobj123 <- FindVariableFeatures(sobj123, selection.method = "vst", nfeatures = 900)
# assume FindVariableFeatures has been run
p1 <- VariableFeaturePlot(sobj123)
top10 <- head(VariableFeatures(sobj123))
p2 <- LabelPoints(plot = p1, points=top10, repel=TRUE)
p2

# # What features account for the most variation within lung samples?
p_dim_loadings <- VizDimLoadings(sobj123, dims = 1:2, reduction = "pca") + ggtitle("Largest sources of variation")
p_dim_loadings
#ggsave(paste0(results_dir, "VizDimLoadings_patients.png"), p_dim_loadings)

DimHeatmap(sobj123, dims = 1:7, cells = 500, balanced = TRUE)


# Check for batch effects?
DimPlot(sobj123, reduction = "pca", group.by = "Sample.Label", alpha=0.5, cols='glasbey')
DimPlot(sobj123, reduction = "pca", group.by = "PanCK.PT", alpha=0.5, cols='glasbey')
dim(sobj123@meta.data) #
DimPlot(sobj123, reduction = "pca", group.by = "fov", alpha=0.5, cols='glasbey')


# What features explain most of the lung variance?
library(wesanderson)
colors <- colorRampPalette(wes_palette("GrandBudapest1"))(33)
colors <- colorRampPalette(c("red", "blue", "green", "yellow","orange", "brown","grey80", "purple"))(30)

colors
colors_tumor <- colorRampPalette(c("red","orange", "peachpuff", "brown","grey80", "purple"))(9)
colors_io <- colorRampPalette(c("blue", "green", "grey50", "yellow"))(25)

#
sobj123 <- FindNeighbors(sobj123, dims = 1:20)
# Computes nearest neighbor graph, SNN in sobj@graphs$SCT_nn, SCT_nn. This is used for clustering.
sobj123 <- FindClusters(sobj123, resolution=0.5, cluster.name='seurat_cluster.0.5', verbose = TRUE) # try different values of resolution for more/less clust.
# lowest resolution (0.1 finds 4 communities, 0.3 finds 9 )
sobj123 <- RunUMAP(sobj123, dims = 1:20, min_dist = 0.5, n_neighbors=20) # from 0.1, 30

DimPlot(sobj123, reduction = "umap", group.by = "seurat_cluster.0.5", alpha=0.5, cols=colors)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_full15", alpha=0.5, cols=colors)
DimPlot(sobj123, reduction = "umap", group.by = "InSituType.level1", alpha=0.5)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_15new1", alpha=0.5, cols=colors)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_Labeled", alpha=0.5, cols=colors)

# myclust$clust[myclust$clust == "i"] <- "A549(adeno)"
# myclust$clust[myclust$clust == "h"] <- "HCC78(adeno)"
# myclust$clust[myclust$clust == "d"] <- "H1975(adeno)"
# myclust$clust[myclust$clust == "o"] <- "SKMES1(SCC)"
# myclust$clust[myclust$clust == "g"] <- "H1703(SCC)"
# myclust$clust[myclust$clust == "f"] <- "HOP62(adeno)"
# myclust$clust[myclust$clust == "j"] <- "HOP92(adeno)"

# Assign Labels
df_meta <- sobj@meta.data
df_meta <- df_meta %>%
  mutate(CellType_Detail_Label = case_when(
    CellType_15new8 == 'i' ~ "A549(adeno)",  # Set 'i' to 'A1'
    CellType_15new8 == 'h' ~ "HCC78(adeno)",
    CellType_15new8 == 'd' ~ "H1975(adeno)",
    CellType_15new8 == 'o' ~ "SKMES1(SCC)",
    CellType_15new8 == 'g' ~ "H1703(SCC)",
    CellType_15new8 == 'b' ~ "H596(adenoSq)",
    CellType_15new8 == 'e' ~ "H596(adenoSq)",
    #CellType_full15 == 'mix' ~ 'B3',
    #CellType_full15 == 'g and f?' ~ 'C1',
    CellType_15new8 == 'f' ~ "HOP62(adeno)", # also some g
    CellType_15new8 == 'j' ~ "HOP92(adeno)",
    CellType_15new8 == 'c' ~ 'Tumor1a',
    CellType_15new8 == 'c_1' ~ 'Tumor1',
    CellType_15new8 == 'c_2' ~ 'Tumor2',
    CellType_15new8 == 'c_3' ~ 'Tumor3',
    CellType_15new8 == 'c_4' ~ 'Tumor4',
    CellType_15new8 == 'c_5' ~ 'Tumor5',
    CellType_15new8 == 'c_6' ~ 'Tumor6',
    CellType_15new8 == 'c_7' ~ 'Tumor7',
    CellType_15new8 == 'c_8' ~ 'Tumor8',
    CellType_15new8 == 'n' ~ 'Tumor2',
    TRUE ~ CellType_15new8)) %>%
  mutate(CellType_Labeled = case_when(
      CellType_15new8 == 'i' ~ 'A1',  # Set 'i' to 'A1'
      CellType_15new8 == 'h' ~ 'A2',
      CellType_15new8 == 'd' ~ 'A3',
      CellType_15new8 == 'o' ~ 'B1',
      CellType_15new8 == 'g' ~ 'B2',
      CellType_15new8 == 'b' ~ "B3",
      CellType_15new8 == 'e' ~ "B3",
      #CellType_full15 == 'mix' ~ 'B3',
      #CellType_full15 == 'g and f?' ~ 'C1',
      CellType_15new8 == 'f' ~ 'C1', # also some g
      CellType_15new8 == 'j' ~ 'C2',
      CellType_15new8 == 'c' ~ 'Tumor1a',
      CellType_15new8 == 'c_1' ~ 'Tumor1',
      CellType_15new8 == 'c_2' ~ 'Tumor2',
      CellType_15new8 == 'c_3' ~ 'Tumor3',
      CellType_15new8 == 'c_4' ~ 'Tumor4',
      CellType_15new8 == 'c_5' ~ 'Tumor5',
      CellType_15new8 == 'c_6' ~ 'Tumor6',
      CellType_15new8 == 'c_7' ~ 'Tumor7',
      CellType_15new8 == 'c_8' ~ 'Tumor8',
      CellType_15new8 == 'n' ~ 'Tumor2',
      TRUE ~ CellType_15new8)
)

df_meta <- df_meta %>%
  mutate(CellType_Labeled.2 = case_when(
    CellType_Labeled == 'A1' ~ 'Tumor Cell Line',  # Set 'i' to 'A1'
    CellType_Labeled == 'A2' ~ 'Tumor Cell Line',
    CellType_Labeled == 'A3' ~ 'Tumor Cell Line',
    CellType_Labeled == 'B1' ~ 'Tumor Cell Line',
    CellType_Labeled == 'B2' ~ 'Tumor Cell Line',
    CellType_Labeled == 'B3' ~ 'Tumor Cell Line',
    CellType_Labeled == 'C1' ~ 'Tumor Cell Line',
    CellType_Labeled == 'C2' ~ 'Tumor Cell Line', # also some g
    CellType_Labeled == 'a' ~ 'Tumor Cell Line',
    CellType_Labeled == 'b' ~ 'Tumor Cell Line',
    CellType_Labeled == 'k' ~ 'Tumor Cell Line',
    CellType_Labeled == 'l' ~ 'Tumor Cell Line',
    CellType_Labeled == 'm' ~ 'Tumor Cell Line',
    CellType_Labeled == 'e' ~ 'Tumor Cell Line',
    CellType_Labeled == 'Tumor1a' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor1' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor2' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor3' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor4' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor5' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor6' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor7' ~ 'Patient Tumor',
    CellType_Labeled == 'Tumor8' ~ 'Patient Tumor',
    TRUE ~ CellType_Labeled
  ))

df_meta <- df_meta %>%
  mutate(CellType_Labeled.1 = case_when(
    CellType_Labeled.2 == 'Tumor Cell Line' ~ 'Tumor',  # Set 'i' to 'A1'
    CellType_Labeled.2 == 'Patient Tumor' ~ 'Tumor',
    CellType_Labeled.2 == 'Endothelial' ~ 'Stroma',
    TRUE ~ 'Stroma'
  ))

sobj <- AddMetaData(sobj, df_meta[,c('cell_id', 'CellType_Labeled', 'CellType_Labeled.1', 'CellType_Labeled.2', 'CellType_Detail_Label')])
sobj123 <- AddMetaData(sobj123, df_meta[c('cell_id', 'CellType_Labeled', 'CellType_Labeled.1',  'CellType_Labeled.2')])


df_meta <- sobj@meta.data
df_meta <- df_meta %>%
  mutate(CellType_Detail_Label5 = case_when(
    CellType_15new5 == 'i' ~ "A549(adeno)",  # Set 'i' to 'A1'
    CellType_15new5 == 'h' ~ "HCC78(adeno)",
    CellType_15new5 == 'd' ~ "H1975(adeno)",
    CellType_15new5 == 'o' ~ "SKMES1(SCC)",
    CellType_15new5 == 'g' ~ "H1703(SCC)",
    CellType_15new5 == 'b' ~ "H596(adenoSq)",
    CellType_15new5 == 'e' ~ "H596(adenoSq)",
    #CellType_full15 == 'mix' ~ 'B3',
    #CellType_full15 == 'g and f?' ~ 'C1',
    CellType_15new5 == 'f' ~ "HOP62(adeno)", # also some g
    CellType_15new5 == 'j' ~ "HOP92(adeno)",
    CellType_15new5 == 'c' ~ 'Tumor1a',
    CellType_15new5 == 'c_1' ~ 'Tumor1',
    CellType_15new5 == 'c_2' ~ 'Tumor2',
    CellType_15new5 == 'c_3' ~ 'Tumor3',
    CellType_15new5 == 'c_4' ~ 'Tumor4',
    CellType_15new5 == 'c_5' ~ 'Tumor5',
    CellType_15new5 == 'c_6' ~ 'Tumor6',
    CellType_15new5 == 'c_7' ~ 'Tumor7',
    CellType_15new5 == 'c_8' ~ 'Tumor8',
    CellType_15new5 == 'n' ~ 'Tumor2',
    TRUE ~ CellType_15new5)) %>%
  mutate(CellType_Labeled5 = case_when(
    CellType_15new5 == 'i' ~ 'A1',  # Set 'i' to 'A1'
    CellType_15new5 == 'h' ~ 'A2',
    CellType_15new5 == 'd' ~ 'A3',
    CellType_15new5 == 'o' ~ 'B1',
    CellType_15new5 == 'g' ~ 'B2',
    CellType_15new5 == 'b' ~ "B3",
    CellType_15new5 == 'e' ~ "B3",
    #CellType_full15 == 'mix' ~ 'B3',
    #CellType_full15 == 'g and f?' ~ 'C1',
    CellType_15new5 == 'f' ~ 'C1', # also some g
    CellType_15new5 == 'j' ~ 'C2',
    CellType_15new5 == 'c' ~ 'Tumor1a',
    CellType_15new5 == 'c_1' ~ 'Tumor1',
    CellType_15new5 == 'c_2' ~ 'Tumor2',
    CellType_15new5 == 'c_3' ~ 'Tumor3',
    CellType_15new5 == 'c_4' ~ 'Tumor4',
    CellType_15new5 == 'c_5' ~ 'Tumor5',
    CellType_15new5 == 'c_6' ~ 'Tumor6',
    CellType_15new5 == 'c_7' ~ 'Tumor7',
    CellType_15new5 == 'c_8' ~ 'Tumor8',
    CellType_15new5 == 'n' ~ 'Tumor2',
    TRUE ~ CellType_15new5)
  )

df_meta <- df_meta %>%
  mutate(CellType_Labeled5.2 = case_when(
    CellType_Labeled5 == 'A1' ~ 'Tumor Cell Line',  # Set 'i' to 'A1'
    CellType_Labeled5 == 'A2' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'A3' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'B1' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'B2' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'B3' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'C1' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'C2' ~ 'Tumor Cell Line', # also some g
    CellType_Labeled5 == 'a' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'b' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'k' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'l' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'm' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'e' ~ 'Tumor Cell Line',
    CellType_Labeled5 == 'Tumor1a' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor1' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor2' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor3' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor4' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor5' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor6' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor7' ~ 'Patient Tumor',
    CellType_Labeled5 == 'Tumor8' ~ 'Patient Tumor',
    TRUE ~ CellType_Labeled5
  ))

df_meta <- df_meta %>%
  mutate(CellType_Labeled5.1 = case_when(
    CellType_Labeled5.2 == 'Tumor Cell Line' ~ 'Tumor',  # Set 'i' to 'A1'
    CellType_Labeled5.2 == 'Patient Tumor' ~ 'Tumor',
    CellType_Labeled5.2 == 'Endothelial' ~ 'Stroma',
    TRUE ~ 'Stroma'
  ))

sobj <- AddMetaData(sobj, df_meta[,c('cell_id', 'CellType_Labeled5', 'CellType_Labeled5.1', 'CellType_Labeled5.2', 'CellType_Detail_Label5')])

sobj123_tumor <- subset(sobj123, subset = CellType_Labeled5.1 %in% c('Tumor'))
sobj12_tumor <- subset(sobj123_tumor, subset = Run_Tissue_name %in% c('PSR-01', 'PSR-02')
                       & CellType_Labeled5.2 == 'Patient Tumor' & PatientID == 6)
colors
Idents(sobj12_tumor) <- "CellType_Detail_Label"
colors_tumor <- colorRampPalette(c("red","orange", "yellow", "green", "blue", "purple"))(8)
DimPlot(sobj12_tumor, reduction = "umap", pt.size = 0.5, group.by = "CellType_Detail_Label", alpha=0.5 , cols=colors_tumor)
DimPlot(sobj12_tumor, reduction = "umap", pt.size = 0.5, group.by = "CellType_Detail_Label5", alpha=0.5, cols=colors_tumor)
ImageDimPlot(sobj12_tumor, fov="PSR01", dark.background = FALSE, group.by = "CellType_Detail_Label5", cols=colors_tumor)
# why doesn't fov PSR02 work on sobj12_tumor
ImageDimPlot(sobj12_tumor, fov="PSR02", dark.background = FALSE, 
             group.by = "CellType_Detail_Label5", cols=colors_tumor)

ImageDimPlot(sobj12_tumor, fov="PSR01", dark.background = FALSE, group.by = "CellType_Detail_Label", cols=colors_tumor)
# why doesn't fov PSR02 work on sobj12_tumor
ImageDimPlot(sobj12_tumor, fov="PSR02", dark.background = FALSE, 
             group.by = "CellType_Detail_Label", cols=colors_tumor)


# todo.
# export table of FindMarkers for each of 5 or 8 tumor clusters. 
# FeaturePlots on umap and on image, in grid. 
# decide which genes; lineage mkrs and also the findmarkers
Idents(sobj12_tumor) <- "CellType_Detail_Label"
sobj12_tumor <- PrepSCTFindMarkers(sobj12_tumor, assay = "SCT", verbose = TRUE)
mkrs8 <- FindAllMarkers(sobj12_tumor, verbose = TRUE)
mkrs8
top_genes <- mkrs8 %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 3) %>%
  ungroup()

top_genes <- top_genes %>%
  arrange(cluster, avg_log2FC)

top_genes$gene

# for each cluster, featureplot the top 3
features_to_plot <- top_genes$gene[top_genes$cluster == "Tumor1"]
features_to_plot
for (clus in unique(top_genes$cluster)) {
  features_to_plot <- top_genes$gene[top_genes$cluster == clus]
  feature_plots <- FeaturePlot(sobj12_tumor, features = features_to_plot, ncol = 3,
                               combine=TRUE)
  feature_plots
  ggsave(paste0(results_dir, "/featureplot_", clus, ".png"), last_plot())
}


################## 5
Idents(sobj12_tumor) <- "CellType_Detail_Label5"
sobj12_tumor <- PrepSCTFindMarkers(sobj12_tumor, assay = "SCT", verbose = TRUE)
mkrs5 <- FindAllMarkers(sobj12_tumor, verbose = TRUE)
mkrs5
top_genes5 <- mkrs5 %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 3) %>%
  ungroup()

top_genes5 <- top_genes5 %>%
  arrange(cluster, avg_log2FC)

top_genes5$gene

features_to_plot <- top_genes5$gene[top_genes5$cluster == "Tumor1"]
features_to_plot
for (clus in unique(top_genes5$cluster)) {
  features_to_plot <- top_genes5$gene[top_genes5$cluster == clus]
  feature_plots <- FeaturePlot(sobj12_tumor, features = features_to_plot, ncol = 3,
                               combine=TRUE)
  #feature_plots
  ggsave(paste0(results_dir, "/featureplot_5clust_", clus, ".png"), last_plot())
}




Idents(s_new) <- "Sample.Label"
DoHeatmap(s_new, features = VariableFeatures(s_new)[1:50])
ggsave(paste0(results_dir, "/heatmap_Pt2_top50_heatmap.png"), last_plot())
FeaturePlot(s_new, features = gentles_adeno, ncol = 2)
RidgePlot(s_new, features = gentles_adeno, ncol = 2)

mkrs <- FindAllMarkers(sobj12_tumor, ident.1 = "Tumor1", ident.2 = "Tumor2", min.pct = 0.25)


features_to_plot <- 
feature_plots <- FeaturePlot(seurat_obj, features = features_to_plot, nrow = 3, ncol = 3)





# exclude cells with very low counts
table(df_meta$CellType_Labeled)
colors36 <- colorRampPalette(c("blue", "green", "yellow","orange", "brown","grey80", "purple", "red"))(36)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_Labeled", alpha=0.5, cols=colors36)

table(df_meta$CellType_Labeled.2)
num_colors <- length(unique(df_meta$CellType_Labeled.2))
num_colors
colors.2 <- colorRampPalette(c("blue", "green", "yellow","orange", "brown","grey80", "purple", "red"))(num_colors)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_Labeled.2", alpha=0.5, cols=colors.2)

num_colors <- length(unique(df_meta$CellType_Labeled.2))
num_colors
colors.1 <- colorRampPalette(c("blue", "green", "yellow","orange", "brown","grey80", "purple", "red"))(num_colors)
DimPlot(sobj123, reduction = "umap", group.by = "CellType_Labeled.1", alpha=0.5, cols=colors.1)

# checkpoint
saveRDS(sobj123, file = paste0(results_dir, "/sobj123_w_InSituType_UMAP.RDS"))
saveRDS(negmeans_vec, file = paste0(results_dir, "/negmeans_vec.RDS"))
saveRDS(merged_mtx, file = paste0(results_dir, "/merged_mtx.RDS"))
negmeans_vec <- readRDS(file = paste0(results_dir, "/negmeans_vec.RDS"))

# resume here
sobj123 <- readRDS(file = paste0(results_dir, "/sobj123_w_InSituType_UMAP.RDS"))
merged_mtx <- readRDS(file = paste0(results_dir, "/merged_mtx.RDS"))

# seruat_cluster.0.2 e.g. is now at attribute in metadata
p_clusters <- DimPlot(sobj123, group.by='seurat_cluster.0.5')
p_clusters
DimPlot(sobj123, group.by='Sample.ID')


# FEATURE Plots
features <- c("PTPRC", "EPCAM", "CD3E")
cell_types_tumor <- c("c", "n", "c_1", "c_2", "c_3", "c_4")

# gen lists of cell ids for highlighting
tumor_c_cells <- colnames(sobj123)[sobj123$CellType_full15 == "c"]
b_cells <- colnames(sobj123)[sobj123$CellType_full15 %in% "B.cell"]
t_cells <- colnames(sobj123)[sobj123$CellType_full15 %in% c("T.cell.CD8", "T.cell.CD4")]
endo_cells <- colnames(sobj123)[sobj123$CellType_full15 %in% "Endothelial"]
macrophage_cells <- colnames(sobj123)[sobj123$CellType_full15 %in% "Macrophage"]

DimPlot2(sobj123, cells.highlight = tumor_c_cells)

#ImageDimPlot(sobj123, fov="PSR01", cells = c(tumor_c_cells, macrophage_cells, endo_cells)) 

# 
# 
# 
# Idents(sobj123) <- "CellType_Labeled"
# ImageDimPlot(sobj123, fov="PSR01") 
# RidgePlot(sobj123, features =  c("PTPRC"), ncol = 3)
# 
# VlnPlot(sobj123, features =  c("PTPRC"))
# 
# FeaturePlot(sobj123, features = c("PTPRC") )

#######################################
# Semi Supervised clustering 
# Obs: 15 additional clusters (on top of Immune cosMx panel) does reasonably well
#    but cluster 'c' is very large and needs to be further subdivided

# what we ran last time (2.5 hrs.)
#######################################
# 
# semisup_full <- run_insitu_type2(sobj,
#                              negmeans_vec,
#                              n_clusts = 10, # number *in addition* to labeled clusters
#                              title_suffix="PatientSamples_10clusts",
#                              group.by="Run_Tissue_name",
#                              insitu_matrix=insitu_matrix_cosMx_IO,
#                              results_dir=results_dir)

semisup_full15 <- run_insitu_type2(sobj,
                                 negmeans_vec,
                                 n_clusts = 15, # number *in addition* to labeled clusters
                                 title_suffix="PatientSamples_15clusts",
                                 group.by="Run_Tissue_name",
                                 insitu_matrix=insitu_matrix_cosMx_IO,
                                 results_dir=results_dir)

sobj <- AddMetaData(sobj, metadata = data.frame(CellType_full15 = semisup_full15$clust))
sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_full15 = semisup_full15$clust))
saveRDS(semisup_full15, file = paste0(results_dir, "/semisup_full15.RDS"))
table(sobj$CellType_full15)

semisup_full15 <- readRDS(file = paste0(results_dir, "/semisup_full15.RDS"))
length(semisup_full15$clust)

#refine the clusters, for tumor and stroma subtypes
counts_matrix <- GetAssayData(object = sobj, assay = "SCT", layer="counts")
counts <- counts_matrix %>%
  as.matrix() %>%
  t()

dim(counts)
cells_to_keep <- rownames(counts)
length(semisup_full15$clust)
semisup_full15$clust[1:4]

semisup_full15$clust <- semisup_full15$clust[names(semisup_full15$clust) %in% cells_to_keep]
semisup_full15$prob <- semisup_full15$prob[names(semisup_full15$prob) %in% cells_to_keep]
#filtered_df <- df[rownames(df) %in% rowname_list, ]
semisup_full15$logliks <- semisup_full15$logliks[rownames(semisup_full15$logliks) %in% cells_to_keep,]
filtered_matrix <- semisup_full15$logliks[rownames(semisup_full15$logliks) %in% cells_to_keep, ]
dim(filtered_matrix)
rownames(filtered_matrix)
semisup_full15$logliks <- filtered_matrix
#length(semisup_full15$clust)
names(semisup_full15$clust)
names(semisup_full15$prob)
head(semisup_full15$logliks)
semisup_full15$clust[1:4]
semisup_full15$logliks[1:4] # NAs is problem.

#newclusts4 <- newclusts
newclusts8 <- refineClusters(
  logliks = semisup_full15$logliks,
  #merges = c("fibroblast" = "stroma", "endothelial" = "stroma"),
  #to_delete = c("a"),
  # subclustering via refineClusters is not recommended for semi-supervised
  # results
  # merge = c("T.cell.CD4"="T-cell",
  #           "T.cell.CD8"="T-cell",
  #           "T.cell.regulatory"= "T-cell",
  #           "Plasmablast" = "Plasma-cell", # to merge, need new name
  #           "Plasma" = "Plasma-cell",
  #           "Plasmacytoid.dendritic.cell" = "pDC", # rename
  #           "Dendritic.cell" = "DC"
  #           ),
  subcluster = list("c"=8),
  counts = counts,
  neg = negmeans_vec
) 
str(newclusts8)
sobj <- AddMetaData(sobj, metadata = data.frame(CellType_15new8 = newclusts8$clust))
sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_15new8 = newclusts8$clust))


#################################
# How many subclusters to make: Try 5
#####
newclusts5 <- refineClusters(
  logliks = semisup_full15$logliks,
  subcluster = list("c"=5),
  counts = counts,
  neg = negmeans_vec
) 
str(newclusts5)

options(max.print = 40)
table(newclusts5$clust)
sobj <- AddMetaData(sobj, metadata = data.frame(CellType_15new5 = newclusts5$clust))
sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_15new5 = newclusts5$clust))




# rename clusters

#newclusts8$clust <- factor(newclusts8$clust, levels = 1:8, 
#                           labels = c("Tumor1", "Tumor2", "Tumor3", "Tumor4", "Stroma1", "Stroma2", "Stroma3", "Stroma4"))




# semisup_full20 <- run_insitu_type2(sobj,
#                                    negmeans_vec,
#                                    n_clusts = 20, # number *in addition* to labeled clusters
#                                    title_suffix="PatientSamples_20clusts",
#                                    group.by="Run_Tissue_name",
#                                    insitu_matrix=insitu_matrix_cosMx_IO,
#                                    results_dir=results_dir)

# get anchors for each cell type
#data.frame(CellType_full20 = semisup_full20$clust)
# sobj <- AddMetaData(sobj, metadata = data.frame(CellType_full20 = semisup_full20$clust))
# sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_full20 = semisup_full20$clust))
# saveRDS(semisup_full20, file = paste0(results_dir, "/semisup_full20.RDS"))

# 
# semisup5 <- run_insitu_type2(sobj123,
#                              negmeans_vec,
#                              n_clusts = 10, # number *in addition* to labeled clusters
#                              title_suffix="PatientSamples_10clusts",
#                              group.by="Run_Tissue_name",
#                              insitu_matrix=insitu_matrix_cosMx_IO,
#                              results_dir=results_dir)
# 
# sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_IO5 = semisup5$clust))
# saveRDS(semisup5, file = paste0(results_dir, "/semisup5.RDS"))
# 
# insitu_cell_typing_plots(sobj, 
#                          insitutype_results=newclusts, 
#                          results_dir, 
#                          celltype_colors,
#                          group.by = "Run_Tissue_name", 
#                          title_suffix = "InSituType_full20"
# )

colors37 <- colorRampPalette(c("darkblue", "darkgreen", "yellow","darkorange", "tan","grey80", "darkcyan", "purple", "darkred"))(37)

# # cellines <-  c("adeno",  "adeno" ,"adeno",  "SCC",  "SCC",   "adenoSquamous",   "adeno",    "adeno")
#names(cellines) <- c("HOP92",  "HOP62" ,"H596",  "H1703",  "SKMES1",   "H1975",   "HCC78",    "A549")
# ideally, here we need named colors that match the clusts names
table(newclusts8$clust)
myclust <- newclusts8 
myclust <- newclusts5

# rename the vectors 
myclust$clust[myclust$clust == "i"] <- "A549"
myclust$clust[myclust$clust == "h"] <- "HCC78"
myclust$clust[myclust$clust == "d"] <- "H1975"
myclust$clust[myclust$clust == "o"] <- "SKMES1"
myclust$clust[myclust$clust == "g"] <- "H1703"
myclust$clust[myclust$clust == "b"] <- "H596" 
myclust$clust[myclust$clust == "e"] <- "H596" 
myclust$clust[myclust$clust == "f"] <- "HOP62"
myclust$clust[myclust$clust == "j"] <- "HOP92"
myclust$clust[myclust$clust == "c"] <- "Tumor1a"
myclust$clust[myclust$clust == "c_1"] <- "Tumor1"
myclust$clust[myclust$clust == "c_2"] <- "Tumor2"
myclust$clust[myclust$clust == "c_3"] <- "Tumor3"
myclust$clust[myclust$clust == "c_4"] <- "Tumor4"
myclust$clust[myclust$clust == "c_5"] <- "Tumor5"
myclust$clust[myclust$clust == "c_6"] <- "Tumor6"
myclust$clust[myclust$clust == "c_7"] <- "Tumor7"
myclust$clust[myclust$clust == "c_8"] <- "Tumor8"

myclust$clust[myclust$clust == "Plasmacytoid.dendritic.cell"] <- "pDC"
myclust$clust[myclust$clust == "Dendritic.cell"] <- "DC"
myclust$clust[myclust$clust == "Plasmablast"] <- "Plasma"
myclust$clust[myclust$clust == "T.cell.regulatory"] <- "T-cell"
myclust$clust[myclust$clust == "T.cell.CD4"] <- "T-cell"
myclust$clust[myclust$clust == "T.cell.CD8"] <- "T-cell"

# apply same changes to myclust$profiles
cols <- colnames(myclust$profiles)
cols <- replace(cols, cols == "i", "A549")
cols <- replace(cols, cols == "h", "HCC78")
cols <- replace(cols, cols == "d", "H1975")
cols <- replace(cols, cols == "o", "SKMES1")
cols <- replace(cols, cols == "g", "H1703")
cols <- replace(cols, cols == "b", "H596")
cols <- replace(cols, cols == "e", "H596")
cols <- replace(cols, cols == "f", "HOP62")
cols <- replace(cols, cols == "j", "HOP92")
cols <- replace(cols, cols == "c", "Tumor1a")
cols <- replace(cols, cols == "c_1", "Tumor1")
cols <- replace(cols, cols == "c_2", "Tumor2")
cols <- replace(cols, cols == "c_3", "Tumor3")
cols <- replace(cols, cols == "c_4", "Tumor4")
cols <- replace(cols, cols == "c_5", "Tumor5")
cols <- replace(cols, cols == "c_6", "Tumor6")
cols <- replace(cols, cols == "c_7", "Tumor7")
cols <- replace(cols, cols == "c_8", "Tumor8")
cols <- replace(cols, cols == "Plasmacytoid.dendritic.cell", "pDC")
cols <- replace(cols, cols == "Dendritic.cell", "DC")
cols <- replace(cols, cols == "Plasmablast", "Plasma")
cols <- replace(cols, cols == "T.cell.regulatory", "T-cell")
cols <- replace(cols, cols == "T.cell.CD4", "T-cell")
cols <- replace(cols, cols == "T.cell.CD8", "T-cell")
colnames(myclust$profiles) <- cols
colnames(myclust$profiles)

tumorclust <- myclust 
tumorclust4 <- tumorclust
tumorclust6 <- tumorclust

# now filter myclust to only cols we want
myclust$profiles <- myclust$profiles[, colnames(myclust$profiles) %in%
                                       c("A549", "HCC78", "H1975", "SKMES1", 
                                         "H1703", "HOP62", "HOP92",
                                          "H596", "H596")]
# merge H596's into one (better to do upstream but oh well)
                                     # skipped a, k, l, m, n    

# patient samples, keep more
tumorclust6$profiles <- tumorclust6$profiles[, colnames(tumorclust$profiles) %in%
                                          c("A549", "HCC78", "H1975", "SKMES1", 
                                            "H1703", "HOP62", "HOP92",
                                            "H596", 
                                            "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6")]

                                     
  #"Tumor1a", "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6", "Tumor7", "Tumor8", "pDC", "DC", "Plasma", "T-cell")]

head(myclust$profiles)
head(tumorclust$profiles)
heatmap(sweep(myclust$profiles, 1, pmax(apply(myclust$profiles, 1, max), .2), "/"), scale="none",
        main = "CellLine Signatures",
)

mtx <- (sweep(myclust$profiles, 1,pmax(apply(myclust$profiles, 1, max), .2)
          , "/"))
# which rows have low totals
threshold <- 0.4
# Remove rows with a row sum less than the threshold
filtered_matrix <- mtx[rowSums(mtx) >= threshold, ]
filtered_matrix <- mtx[apply(mtx, 1, max) >= threshold, ]
dim(filtered_matrix)

heatmap(filtered_matrix, scale="none",
        main = "CellLine Signatures",
)

### Reduce heatmap
# Patient 6
tumorclust6$profiles <- tumorclust6$profiles[, colnames(tumorclust6$profiles) %in%
                                             c("A549", "HCC78", "H1975", "SKMES1", 
                                               "H1703", "HOP62", "HOP92",
                                               "H596", 
                                               "Tumor1", "Tumor3", "Tumor4", "Tumor6")]
# rename: 
cols <- colnames(tumorclust6$profiles)
cols
cols <- replace(cols, cols == "Tumor1", "Pt6.post")
cols <- replace(cols, cols == "Tumor3", "Pt6.pre.Reg2")
cols <- replace(cols, cols == "Tumor4", "Pt6.pre.Reg1")
cols <- replace(cols, cols == "Tumor6", "Pt6.pre.Reg3")
colnames(tumorclust6$profiles) <- cols


mtx <- (sweep(tumorclust6$profiles, 1,pmax(apply(tumorclust6$profiles, 1, max), .2) , "/"))
# which rows have low totals
threshold <- 0.3
# Remove rows with a row sum less than the threshold
filtered_matrix <- mtx[rowSums(mtx) >= threshold, ]
filtered_matrix <- mtx[apply(mtx, 1, max) >= threshold, ]
dim(filtered_matrix)

# Set appropriate margins to avoid cutting off axis labels
par(mar = c(1, 1, 1, 1))
heatmap(filtered_matrix, scale="none",
        main = "Patient 6 CellLine Signatures",
        margins = c(5, 5),
        
)

# Patient 6
tumorclust4$profiles <- tumorclust4$profiles[, colnames(tumorclust4$profiles) %in%
                                               c("A549", "HCC78", "H1975", "SKMES1", 
                                                 "H1703", "HOP62", "HOP92",
                                                 "H596", 
                                                 "Tumor2", "Tumor5")]
# rename: 
cols <- colnames(tumorclust4$profiles)
cols
cols <- replace(cols, cols == "Tumor2", "Pt4.post")
cols <- replace(cols, cols == "Tumor5", "Pt4.pre")
colnames(tumorclust4$profiles) <- cols


mtx <- (sweep(tumorclust4$profiles, 1,pmax(apply(tumorclust4$profiles, 1, max), .2) , "/"))
# which rows have low totals
threshold <- 0.3
# Remove rows with a row sum less than the threshold
filtered_matrix <- mtx[rowSums(mtx) >= threshold, ]
filtered_matrix <- mtx[apply(mtx, 1, max) >= threshold, ]
dim(filtered_matrix)

# Set appropriate margins to avoid cutting off axis labels
par(mar = c(1, 1, 1, 1))
heatmap(filtered_matrix, scale="none",
        main = "Patient 4 CellLine Signatures",
        margins = c(5, 5),
)



SeuratExtend::Heatmap(filtered_matrix,
                       scale = "none", main = "CellLine Signatures", 
                      filename = paste0(results_dir, "/Heatmap_CellLines.png"))
                      
#ggsave(paste0(results_dir, "/Heatmap_CellLines.png"), width = 10, height = 10, dpi = 300)

sobj <- AddMetaData(sobj, metadata = data.frame(CellType_Labeled.3 = myclust$clust))

df_meta <- sobj@meta.data
df_meta<- df_meta %>%
  mutate(Sample.Label.4 = case_when(
    Diagnosis == "" ~ df_meta$CellType_Labeled,
    TRUE ~ Diagnosis
  ))
table(df_meta$Sample.Label.4)
sobj <- AddMetaData(sobj, metadata = data.frame(Sample.Label.4 = df_meta$Sample.Label.4))

# specific color assingments:
colors_tumor1 <- c("#FF0000", "#FF3300", "#FF6700", 
                   "#A020F0", "#BB8BD9", "brown",
                   "darkorange", "#FF9A00",      # 
                   "#FFD3A1", "#EEB99E", "#D28271", "#B54B44", "#FFC268", 
                   "#FFB22E", "#B67070" , "#FFC268", "#FFB22E")
names(colors_tumor1) <- c("A1","A2","A3","B1","B2","B3","C1","C2", "Tumor1a",
                          "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6", "Tumor7", 
                          "Tumor8")

colors_tumor_cpa <- c("#FF0000", "#FF3300", "#FF6700", 
                   "#A020F0", "#BB8BD9", "brown",
                   "darkorange", "#FF9A00",      # 
                   "#FFD3A1", "#EEB99E", "#D28271", "#B54B44", "#FFC268", 
                   "#FFB22E", "#B67070" , "#FFC268", "#FFB22E")
names(colors_tumor_cpa) <- c("A1","A2","A3","B1","B2","B3","C1","C2", "Tumor1a",
                          "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6", "Tumor7", 
                          "Tumor8")

iocolors
# but change some colors
#iocolors["neutrophil"] <- "#9966CC" 
#iocolors["plasmaplast"] <- "#3399CC"

# 1. seqalong to make sure all items get a color
cols <- colors37[seq_along(unique(sobj@meta.data$CellType_Labeled))]
names(cols) <- unique(sobj@meta.data$CellType_Labeled)
cols
# 2. match up named colors in iocolors + colors_tumor1
colors_named <- c(iocolors, colors_tumor1)
cols[is.element(names(cols), names(colors_named))] <- colors_named[names(cols)[is.element(names(cols), names(colors_named))]]

# check all celltypes are in the colorset
length(unique(names(cols)))
length(unique(sobj@meta.data$CellType_Labeled))
setdiff(unique(newclusts8$clust), unique(names(cols)))
setdiff(unique(names(cols)), unique(newclusts8$clust))

table(sobj@meta.data$CellType_Labeled)
table(newclusts8$clust)

#### Flightpath


# Save 


#colors37 <- c("grey", "red", "blue", "green", "purple", "orange", "yellow", "black", "brown", "pink", "cyan", "magenta", "darkgreen", "darkblue", "darkred", "darkorange", "darkyellow", "darkgrey", "darkpurple", "darkcyan")
insitu_cell_typing_plots(sobj, 
                         insitutype_results=newclusts8, 
                         results_dir, 
                         colorset = cols,
                         group.by = "Run_Tissue_name", 
                         title_suffix = "CellTypes_High_Detail",
                         "CellType_Labeled"
)




##########################
# FLIGHTPATH  ONLY


flightpath <- InSituType::flightpath_layout(logliks = myclust$logliks, 
                                            profiles = myclust$profiles)

# Save the plot as a PNG file
png(paste0(results_dir, paste0("/InSitu_FlightPath_for_celltypes.png")), width = 600, height = 800)  # Open the PNG device
par(mar = c(0,0,0,0))
plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols)
text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
dev.off()


cols = colors37
# define cluster colors:
cols2 <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )
cols <- cols[seq_along(unique(newclusts8$clust))]
names(cols) <- unique(newclusts8$clust)
cols
iocolors
col
flightpath <- InSituType::flightpath_layout(logliks = newclusts8$logliks, 
                                            profiles = newclusts8$profiles)
par(mar = c(0,0,0,0))
plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[newclusts8$clust])
text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
dev.off()

common_names <- intersect(names(semisup_full15$clust), tumor_anchor_cells)
subset_clust <- semisup_full15$clust[common_names]
table(subset_clust)
# mostly "c" and "n", but the 20 more clusters version moved many to NK.cells (wrong)

# what are my annotated tumor cells annotated as?
common_names <- intersect(names(newclusts8$clust), tumor_anchor_cells)
subset_clust <- newclusts8$clust[common_names]
table(subset_clust)
# mostly 'c' and 'n'

# flightpath

cols <- InSituType::colorCellTypes(freqs = table(semisup10$clust), palette = "brewers")
str(cols)
#>  Named chr [1:21] "#8DD3C7" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" ...
#>  - attr(*, "names")= chr [1:21] "f" "b" "a" "c" ...
#>  
colors[c("Endothelial", "Fibroblast", "B.cell")]

# setNames assings color names to a vector of colors

flightpath <- InSituType::flightpath_layout(logliks = semisup10$logliks, profiles = semisup10$profiles)
par(mar = c(0,0,0,0))
plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[semisup10$clust])
text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)


# get Umap
# plot(um, pch = 16, cex = 0.1, col = cols[semisup10$clust])
# for (cell in unique(semisup10$clust)) {
#   text(median(um[semisup10$clust == cell, 1]), median(um[semisup10$clust == cell, 2]), cell, cex = 0.7)  
# }

# spatial
# xy space, by slide or sample
# how to grou semisup10 by sample
sample.groups <- setNames(sobj@meta.data$Sample.ID, sobj@meta.data$cell_id)
sample.groups[1:5]
head(semisup10$clust)
combined_df <- data.frame(semisup10$clust, sample.groups, row.names = names(semisup10$clust))
head(combined_df)

# for each sample.group in combined_df
#   plot the xy points for that sample.group
#   color the points by semisup10$clust

# Assuming sample.groups is a vector and combined_df is your dataframe
unique_groups <- unique(sample.groups)

df_meta <- sobj@meta.data
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
  
  p <- ggplot(df_meta_grp,
              aes(CenterX_global_px, y=CenterY_global_px, color=CellType_IOPlus)) +
    geom_point() +
    scale_color_manual(values = cols[df_meta$CellType_IOPlus]) +
    ggtitle(paste("Sample", group))
  ggsave(paste0(results_dir, "/sample_", group, "_xy.png"), p)
}





### HERE!

# par(mar = c(1,1,1,1))
# p <- plot(df_meta_grp, aes(x=df_meta_grp$CenterX_global_px, y=df_meta_grp$CenterY_global_px))
#      #pch = 16, cex = 0.1, col = cols[semisup10$clust], asp = 1) 
# print(p)
# how to make a grid of plots
#library(patchwork)
# Create a list of ggplot objects
plots <- lapply(1:6, function(i) {
  ggplot(data.frame(x = 1:10, y = rnorm(10, i, 1)), aes(x = x, y = y)) +
    geom_point() +
    ggtitle(paste("Plot", i))
})

# Combine the plots
plot_grid <- plots[[1]] + plots[[2]] + plots[[3]] +
  plots[[4]] + plots[[5]] + plots[[6]] +
  plot_layout(ncol = 3, nrow = 2)

# Print the plot grid
plot_grid

# 




semisupIO_10b <- run_insitu_type2(sobj123,
                       negmeans_vec,
                       n_clust = 3, # IO has 15 clusters, n_clust is how many *new* clusters to make
                       title_suffix="Patient Samples 10 clust",
                       group.by="Run_Tissue_name",
                       insitu_matrix=as.matrix(merged_mtx),
                       results_dir=results_dir)
# Negative control means are all zero. Please check the negative control names.
# Why did this place so many cells in unassigned clusters?
# differences from the really nice one are: mixed matrix (plus Tumor), and CPA included.
# CPA meant most of the 10 were occupied by CPAs. Try
# 1. less clusters
# 2. revert to mixed matrix. 
# 3. [re-add CPA]
# 
# sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_IO = semisupIO_10b$clust))
# insitu_cell_typing_plots(sobj123, 
#                          insitutype_results=semisupIO_10b, 
#                          results_dir, 
#                          celltype_colors,
#                          group.by = "Run_Tissue_name", 
#                          title_suffix = "InSituType_10b"
#                          )

# flightpath

cols <- InSituType::colorCellTypes(freqs = table(semisup10$clust), palette = "brewers")
str(cols)
#>  Named chr [1:21] "#8DD3C7" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" ...
#>  - attr(*, "names")= chr [1:21] "f" "b" "a" "c" ...
cols[c("Endothelial", "Fibroblast", "B.cell")]

flightpath <- InSituType::flightpath_layout(logliks = semisup10$logliks, profiles = semisup10$profiles)
par(mar = c(0,0,0,0))
plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[semisup10$clust])
text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)



plot(um, pch = 16, cex = 0.1, col = cols[semisup10$clust])
for (cell in unique(semisup10$clust)) {
  text(median(um[semisup10$clust == cell, 1]), median(um[semisup10$clust == cell, 2]), cell, cex = 0.7)  
}

# spatial
# xy space, by slide or sample
# how to grou semisup10 by sample
sample.groups <- setNames(sobj@meta.data$Sample.ID, sobj@meta.data$cell_id)
sample.groups[1:5]
head(semisup10$clust)
combined_df <- data.frame(semisup10$clust, sample.groups, row.names = names(semisup10$clust))
head(combined_df)

# for each sample.group in combined_df
#   plot the xy points for that sample.group
#   color the points by semisup10$clust

# Assuming sample.groups is a vector and combined_df is your dataframe
unique_groups <- unique(sample.groups)

df_meta <- sobj@meta.data
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
  
  p <- ggplot(df_meta_grp,
              aes(CenterX_global_px, y=CenterY_global_px, color=CellType_IOPlus)) +
    geom_point() +
    scale_color_manual(values = cols[df_meta$CellType_IOPlus]) +
    ggtitle(paste("Sample", group))
  ggsave(paste0(results_dir, "/sample_", group, "_xy.png"), p)
}




df_meta <- sobj@meta.data
df_meta$cell_type <- semisup10$clust

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

# this function expects 'cell_type' col 
p <- plot_celltypes(df_celltype, "Run_Tissue_name", "")

# compare
class(insitu_matrix)
class(merged_mtx)

# sup <- insitutypeML(x = counts,
#                   neg = negmean,
#                   cohort = NULL, # cohort_data,
#                   #n_clusts = 0, # start generous
#                   #n_anchor_cells = 1000,
#                   #insufficient_anchors_thresh = 20,
#                   reference_profiles = merged_mtx)  

# automatically selecting anchor cells with the best fits to fixed profiles
# The following cell types had too few anchors and so are being removed from consideration: X6_5, X6_3, X3_6, X3_5, X2_4, X2_3, X6_5.1, X6_3.1, X7_4, X4_1, X5_2, X4_6, B.cell, Dendritic.cell, Endothelial, Fibroblast, Macrophage, Mast.cell, Monocyte, Neutrophil, NK.cell, Plasma, Plasmablast, Plasmacytoid.dendritic.cell, T.cell.CD4, T.cell.CD8, T.cell.regulatory
# Error in updateReferenceProfiles(reference_profiles, counts = x, neg = neg,  : 
#                                    No anchors were selected. The algorithm can't run under these conditions. 
#          Solutions include: 1. make anchor selection more generous. 2. select anchors by hand.
# In addition: There were 50 or more warnings (use warnings() to see the first 50)

# 
# heatmap(sweep(semisup10$profiles, 1, pmax(apply(semisup10$profiles, 1, max), .2), "/"), scale = "none",
#         main = "Mean cell type expression profiles")
# 
# # make the flightpath plot
# fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = semisup10, col = cols[semisup10$clust])
# print(fp)
# 
# # refine the clusters, for tumor and stroma subtypes
# newclusts <- refineClusters(
#   logliks = semisup$logliks,
#   #merges = c("fibroblast" = "stroma", "endothelial" = "stroma"),
#   merges = c("macrophages" = "myeloid",  # merge 3 clusters
#              "monocytes" = "myeloid",
#              "mDC" = "myeloid",
#              "T.cell.CD4" = "T.cell",
#              "T.cell.CD8" = "T.cell",
#              
#              "B-cells" = "lymphoid",
#              "a" = "A1",
#              "b" = "A3",
#              "c" = "A1",
#              "d" = "A1",
#              "e" = "Tumor",
#              "f" = "B1",
#              "g" = "A1",
#              "h" = "A1",
#              "i" = "A2/B2/B3",
#              "j" = "A1",
#              ) ,
#   #to_delete = c("a"),
#   # subclustering via refineClusters is not recommended for semi-supervised
#   # results
#   subcluster = list("e"=4),
#   counts = counts,
#   neg = negmean
# ) 
# str(newclusts)
# 

# 
# 
# sup <- insitutype(x = counts,
#                     neg = negmean,
#                     cohort = NULL, # cohort_data,
#                     n_clusts = 0, # dictates supervised
#                     reference_profiles = merged_mtx)   
# 

# ImmuneTumor_safeTME_profilematrix.csv has 693 genes not matching  
# try Lung_HCA_profilematrix

# TODO: provide alias table and rename reference 
# find 'VEGF' in ioprofiles
#ioprofiles['VEGF',]
# 
# sup$clus[1:10]
# length(sup$clus)
# 
# # nice plots
# heatmap(sweep(sup$profiles, 1, pmax(apply(sup$profiles, 1, max), .2), "/"), scale="none",
#         main = "InSituType Clustering - Anchor Tumor Cells",
# )
# 

########################
# HEATMAP
#########################
# CPA and also Patient Tumor Profiles


# Tumor cells only
head(newclusts8$profiles)
head(myclust$profiles)
cols <- colnames(newclusts8$profiles)
class(cols)
cols <- replace(cols, cols == "a", "x")
cols <- replace(cols, cols == "b", "y")

heatmap(sweep(semisup_full15$profiles, 1, pmax(apply(semisup_full15$profiles, 1, max), .2), "/"), scale="none",
          main = "Tumor Signatures",
  )
heatmap(sweep(newclusts8$profiles, 1, pmax(apply(newclusts8$profiles, 1, max), .2), "/"), scale="none",
        main = "Tumor Signatures",
)



heatmap(sweep(semisup_full15$profiles, 1, pmax(apply(semisup_full15$profiles, 1, max), .2), "/"), scale="none",
        main = "Tumor Signatures",
)
# heatmap isn't a ggplot so this doesn't work; Save the plot as a file instead.
ggsave(paste0(results_dir,"/ImageDimPlot_NSCLC_heatmap_allslides.png"), plot = p)


p <- pheatmap::pheatmap(sweep(refprofiles, 1, pmax(apply(refprofiles, 1, max), 0.2), "/"), 
                        col = colorRampPalette(c("white", "darkblue"))(100))

p 
ggsave(paste0(results_dir, "/cosMxIO_refprofiles_heatmap.png"), plot = p, width = 10, height = 10, dpi = 300)


# instead of semisup, use saved metadata
table(sobj123@meta.data$CellType_Labeled)
colnames(sobj123@meta.data)
table(sobj123@meta.data$CellType_Labeled.2)
# want Tumor Cell Line and Patient Tumor
table(sobj123@meta.data$CellType_Labeled.1) # Tumor


# use the saved metadata in Sample.Lable.4
# cellines <-  c("adeno",  "adeno" ,"adeno",  "SCC",  "SCC",   "adenoSquamous",   "adeno",    "adeno")
#names(cellines) <- c("HOP92",  "HOP62" ,"H596",  "H1703",  "SKMES1",   "H1975",   "HCC78",    "A549")
df_meta <- sobj@meta.data
df_meta <- df_meta %>%
  mutate(CellType_Labeled.4 = case_when(
    CellType_Labeled %in% c("A1","A2","A3", "C1","C2") ~ "CellLine Adeno",
    CellType_Labeled %in% c("B1","B2") ~ "CellLine Squamous",
    CellType_Labeled %in% c("B3") ~ "CellLine AdenoSquamous",
    CellType_Labeled %in% c("Plasmacytoid.dendritic.cell") ~ "Dendritic",
    CellType_Labeled %in% c("a", "b", "e","k","l","m") ~ "Unlabeled Pt Tumor",
    CellType_Labeled %in% c("T.cell.CD4", "T.cell.CD8", "T.cell.regulatory") ~ "T-cell",
    TRUE ~ CellType_Labeled ))
table(df_meta$CellType_Labeled.4)
head(df_meta)

sobj <- AddMetaData(sobj, metadata = data.frame(CellType_Label.4 = df_meta$SCellType_Label.4))
table(sobj@meta.data$CellType_Labeled.4)

sobj_tumor <- subset(sobj, subset = CellType_Labeled.1 == "Tumor")

sobj_cpa <- subset(sobj_tumor, subset = Run_Tissue_name == "PSR-CPA")
sobj_pt_tumor <- subset(sobj_tumor, subset = Run_Tissue_name != "PSR-CPA")

colors_tumor2 <- cols2[1:length(colors_tumor1)]
names(colors_tumor2) <- names(colors_tumor1)
colors_tumor2

# colors for CellType_Labeled.4
colors_tumor4 <- cols2
names(colors_tumor4) <- unique(sobj_tumor@meta.data$CellType_Labeled.4)
table(sobj_tumor@meta.data$CellType_Labeled.4)
colors_tumor4
colors_tumor4["Tumor3"] <- "#FF0000"
colors_tumor4["CellLine Squamous"] <- "darkblue"
colors_tumor4["CellLine Adeno"] <- "#940E6A" 
colors_tumor4["CellLine Adena"] <- "#940E6A" 


# pick 

e                           k                           b                          B2 
"#00008B"                   "#00166C"                   "#002C4D"                   "#FFB22E" 
l                          A1                          A3                          C2 
"#00580F"                   "#FF0000"                   "#FF6700"                   "#EEB99E" 
B1                      Tumor3                          C1                          A2 
"#FF9A00"                   "#B67070"                   "#FFD3A1"                   "#FF3300" 
Tumor4                      Tumor1                  T.cell.CD8                      Tumor2 
"#C2A3A3"                   "#B54B44"                   "#FA900F"                   "#A93E3E" 
Neutrophil Plasmacytoid.dendritic.cell                      Tumor6                     NK.cell 
"#9966CC"                   "#DCAB6C"                   "#BB8BD9"                      "pink" 
Tumor5                   Mast.cell                  Fibroblast                  T.cell.CD4 
"#C9C1CE"               "springgreen"                   "#999999"                   "#B5C4C4" 
Dendritic.cell                      Tumor7                      Tumor8                  Macrophage 
"#88B6B6"                   "#AD55E4"                   "#A020F0"                   "#006600" 
Endothelial           T.cell.regulatory                    Monocyte                           a 
"#996633"                   "#475BB7"                   "#33CC00"                   "#8E2BE4" 
B.cell                 Plasmablast                           m                      Plasma 
"darkblue"                   "#3399CC"                   "#940E6A"                   "#8F0735" 


#colors_tumor2 <- colors_tumor1
Idents(sobj_tumor) <- "CellType_Labeled.4"
ImageDimPlot(sobj_tumor, fov="PSR02", group.by = "CellType_Labeled.4", cols="glasbey", 
             size=0.8,
             dark.background=FALSE)
ggsave(paste0(results_dir,"/ImageDimPlot_TumorHeterogeneity.png"), width = 10, height = 10, dpi = 300)

Idents(sobj_tumor) <- "CellType_15new5"
ImageDimPlot(sobj_tumor, fov="PSR02", group.by = "CellType_15new5", cols="glasbey", 
             size=0.8,
             dark.background=FALSE)
ggsave(paste0(results_dir,"/ImageDimPlot_TumorHeterogeneity.png"), width = 10, height = 10, dpi = 300)

table(sobj_tumor@meta.data$CellType_15new8)



# CPA needs a calculated column from CellType_15new8
df_meta <- sobj_cpa@meta.data
df_meta <- df_meta %>%
  mutate(CellType_15new8_Label = case_when(
    CellType_15new8_ %in% c("j") ~ "A(adeno)",
    CellType_15new8 %in% c("i") ~ "B(squamous)",
    CellType_Labeled.4 %in% c("CellLine Squamous") ~ "CellLine Squamous",
    CellType_Labeled.4 %in% c("CellLine AdenoSquamous") ~ "CellLine AdenoSquamous",
    CellType_Labeled.4 %in% c("Dendritic") ~ "Dendritic",
    CellType_Labeled.4 %in% c("Unlabeled Pt Tumor") ~ "Unlabeled Pt Tumor",
    CellType_Labeled.4 %in% c("T-cell") ~ "T-cell",
    TRUE ~ CellType_Labeled.4 ))

Idents(sobj_cpa) <- "CellType_15new8"
ImageDimPlot(sobj_tumor, fov="PSRCPA", group.by = "CellType_15new8",  # cols="colors_tumor1",
             size=0.8,
             dark.background=FALSE)

ImageDimPlot(sobj_tumor, fov="PSR02", group.by = "CellType_15new8",  # cols="colors_tumor1",
             size=0.8,
             dark.background=FALSE)


ImageDimPlot(sobj_tumor, fov="PSR02", group.by = "CellType_15new5",  # cols="colors_tumor1",
             size=0.8,
             dark.background=FALSE)

# pheatmap(sweep(sobj_cpa@assays$RNA@counts, 1, pmax(apply(sobj_cpa@assays$RNA@counts, 1, max), .2), "/"), scale="none",
#         main = "InSituType Clustering - CPA Tumors",
# )


# First, create a matrix using the CalcStats function.
Idents(sobj_tumor) <- "Sample.Label"
genes <- VariableFeatures(sobj_tumor) # all of them
topplot <- CalcStats(sobj_tumor, features = genes, method = "zscore", order = "p", n = 5)
# Generate a basic heatmap.
hm1 <- Heatmap(topplot, lab_fill = "zscore")
hm1
ggsave(paste0(results_dir,"/Heatmap_Tumor_all1.png"), plot = hm1, width = 10, height = 40, dpi = 300)

# zoom in: only CPA 
Idents(sobj_cpa) <- "SPORE"
genes <- VariableFeatures(sobj_cpa)
topplot <- CalcStats(sobj_cpa, features = genes, method = "zscore", order = "p", n = 15)
head(topplot)

# how to map samples to type
#cellines <- c("HOP92",  "HOP62" ,"H596",  "H1703",  "SKMES1",   "H1975",   "HCC78",    "A549")
#names(cellines) <- c("adeno",  "adeno" ,"adeno",  "SCC",  "SCC",   "adenoSquamous",   "adeno",    "adeno")

# in reverse
cellines <-  c("adeno",  "adeno" ,"adeno",  "SCC",  "SCC",   "adenoSquamous",   "adeno",    "adeno")
names(cellines) <- c("HOP92",  "HOP62" ,"H596",  "H1703",  "SKMES1",   "H1975",   "HCC78",    "A549")

# Generate a basic heatmap.
hm1 <- Heatmap(topplot, lab_fill = "zscore", facet_col = cellines) + ggtitle("Cell Pellet Array Signatures by type")
hm1 
ggsave(paste0(results_dir,"/Heatmap_TumorType_cpa.png"), plot = hm1, width = 10, height = 40, dpi = 300)

Idents(sobj_cpa) <- "Diagnosis"
genes <- VariableFeatures(sobj_cpa)
topplot <- CalcStats(sobj_cpa, features = genes, method = "zscore", order = "p", n = 15)
head(topplot)
hm1 <- Heatmap(topplot, lab_fill = "zscore") + ggtitle("NSCLC Type Signatures")
hm1 
ggsave(paste0(results_dir,"/Heatmap_TumorType_cpa_general.png"), plot = hm1, width = 10, height = 20, dpi = 300)

# What changed by patient
Idents(sobj_pt_tumor) <- "PatientID"
sobj4 <- subset(sobj_tumor, subset = PatientID %in% c(4))
Idents(sobj4) <- "Sample.Label2"
genes <- VariableFeatures(sobj4)
topplot <- CalcStats(sobj4, features = genes, method = "zscore", order = "p", n = 15)
head(topplot)
hm1 <- Heatmap(topplot, lab_fill = "zscore") + ggtitle("Patient 4 Differences")
hm1
ggsave(paste0(results_dir,"/Heatmap_Patient4_Differentiators.png"), plot = hm1, width = 10, height = 20, dpi = 300)

# What changed by patient
Idents(sobj_pt_tumor) <- "PatientID"
sobj6 <- subset(sobj_tumor, subset = PatientID %in% c(6))
Idents(sobj6) <- "Sample.Label2"
genes <- VariableFeatures(sobj6)
topplot <- CalcStats(sobj6, features = genes, method = "zscore", order = "p", n = 15)
head(topplot)
hm1 <- Heatmap(topplot, lab_fill = "zscore") + ggtitle("Patient 6 Differences")
hm1
ggsave(paste0(results_dir,"/Heatmap_Patient6_Differentiators.png"), plot = hm1, width = 10, height = 10, dpi = 300)

# What changed by patient 3
Idents(sobj_pt_tumor) <- "PatientID"
sobj3 <- subset(sobj_tumor, subset = PatientID %in% c(3))
Idents(sobj3) <- "Sample.Label2"
genes <- VariableFeatures(sobj3)
topplot <- CalcStats(sobj3, features = genes, method = "zscore", order = "p", n = 15)
head(topplot)
hm1 <- Heatmap(topplot, lab_fill = "zscore") + ggtitle("Patient 3 Differences")
hm1
ggsave(paste0(results_dir,"/Heatmap_Patient3_Differentiators.png"), plot = hm1, width = 10, height = 10, dpi = 300)

# as waterfall

# make sure group by is a factor and ordered. 
# narrow to paired samples for now
sobj_tumor_pairs <- subset(sobj_tumor, subset = (PatientID %in% c(3,4,6) & Sample.ID != "ORG_4_1")
                           | Run_Tissue_name == "PSR-CPA")
table(sobj_tumor_pairs@meta.data$Sample.Label2)
# # calculate new meta 
# df_meta <- sobj_tumor_pairs@meta.data
# df_meta <- df_meta %>%
#   mutate(Sample.Label3 = case_when(Diagnosis != "" ~ paste(Diagnosis, "cell line"), 
#                                    TRUE ~ Sample.Label2))
# df_meta$Sample.Label3 <- factor(df_meta$Sample.Label3, 
#                            levels = c("220904-liver-pre ( 4_6 )",
#                                       "220904-liver-post ( 4_1 )",
#                                       "220901-brain-pre ( 6_3 )",
#                                       "220901-lymph node-post ( 6_5 )",
#                                       "Adenocarcinoma cell line",
#                                       "AdenoSquamous cell line",
#                                       "SCC cell line"))

#table(df_meta$Sample.Label3)
#sobj_tumor_pairs <- AddMetaData(sobj_tumor_pairs, df_meta["Sample.Label3"])
  
  #factor(Sample.Label2, levels = c("CPA", "CPA2", "CPA3", "CPA4", "CPA5", "CPA6", "CPA7", "CPA8", "CPA9", "CPA10", "CPA11", "CPA12", "CPA13", "CPA14", "CPA15", "CPA16", "CPA17", "CPA18", "CPA19", "CPA20", "CPA21", "CPA22", "CPA23", "CPA24", "CPA25", "CPA26", "CPA27", "CPA28", "CPA29", "CPA30", "CPA31", "CPA32", "CPA33", "CPA34", "CPA35", "CPA36", "CPA37", "CPA38", "CPA39", "CPA40", "CPA41", "CPA42", "CPA43", "CPA44", "CPA45", "CPA46", "CPA47", "CPA48", "CPA49", "CPA50", "CPA51", "CPA52", "CPA53", "CPA54", "CPA55", "CPA56", "CPA57", "CPA58", "CPA59", "CPA60", "CPA61", "CPA62", "CPA63", "CPA64", "CPA65", "CPA66", "CPA67", "CPA68", "CPA69", "CPA70", "CPA71", "CPA72", "CPA73", "CPA74", "CPA75", "CPA76", "CPA77", "CPA78", "CPA79", "CPA80", "CPA81", "CPA82", "CPA83", "CPA84", "CPA85", "CPA86", "CPA87", "CPA88", "CPA89", "CPA90", "CPA91", "CPA92", "CPA93", "CPA94", "CPA95", "CPA96", "CPA97

Idents(sobj_tumor_pairs) <- "Sample.Label3"
VlnPlot2(sobj_tumor_pairs, features = mor_genes, pt_size = 0.1, ncol = 1,  group.by = "Sample.Label3") 


# todo: add GO terms

# only Patient 4 and similar CPA 
ggsave(paste0(results_dir,"/Heatmap_TumorType_cpa.png"), plot = hm1, width = 10, height = 40, dpi = 300)

# break into Samle.Labels but group by Diagnosis



# patient 4 and most similar CPA
sobj_pt4 <- subset(sobj_pt_tumor, subset = PatientID %in% c(4) )
Idents(sobj_pt_tumor) <- "Sample.Label2"


# TODO: use patient IDs instead of sample.ID in Sample.Label
# how to narrow in on key findings e.g. only Pt4 and CPA similar
# add gene annotations
# which cell lines are what cancer type


#HERE!

# Example matrix of gene expression data (rows = genes, columns = samples)
#expression_matrix <- GetAssayData(sobj_tumor_pairs, layer="data")
#rownames(expression_matrix) <- paste0("Gene", 1:10)
# colnames(expression_matrix) <- c("Sample_Pre_1", "Sample_Post_1", "Sample_Pre_2", 
#                                  "Sample_Post_2", "Control_1", "Control_2")
# 
# # Sample annotations
# annotation_col <- data.frame(
#   Condition = factor(c("Pre", "Post", "Pre", "Post", "Control", "Control")),
#   row.names = colnames(expression_matrix)
# )

# Create the heatmap with clustering on genes and annotations for samples
pheatmap(expression_matrix,
         annotation_col = annotation_col,         # Annotate the columns (samples)
         cluster_rows = TRUE,                     # Cluster the genes (rows)
         cluster_cols = FALSE,                    # Do not cluster samples (columns)
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Divergent color palette
         fontsize_row = 8, fontsize_col = 10)

# Modify the color theme to range from white to dark green.
Heatmap(toplot, lab_fill = "zscore", color_scheme = c("white", muted("green")))

# 
# iocolors # list of hex codes
# cols <- iocolors[seq_along(unique(sup$clust))]
# names(cols) <- unique(sup$clust)

# 
# fp <- flightpath_plot(flightpath_result = NULL,
#                       insitutype_result = sup,
#                       col = cols[sup$clust])
# print(fp)
# ggsave(paste0(results_dir,"ImageDimPlot_NSCLC_FlightPath_all.png"), plot = fp)
# 
# length(sup$clust)
# sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_ns = sup$clust))
# 
# Idents(sobj) <- "CellType_ns"
# #Idents(sobj) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
# p <- ImageDimPlot(sobj123, fov="PSR01", cols="polychrome", dark.background = FALSE)
# show(p)
# #ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_PSR01.png"), plot = p)
# p <- ImageDimPlot(sobj123, fov="PSR02", cols="polychrome", dark.background = FALSE)
# #ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_PSR02.png"), plot = p)
# p <-ImageDimPlot(sobj123, fov="PSROrganoids", cols="polychrome", dark.background = FALSE)
# #ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_Organoids.png"), plot = p)
# 
# 

#####################################
# FindAllMarkers
#####################################
#sobj <- PrepSCTFindMarkers(sobj, assay = "SCT", verbose = TRUE)
Idents(sobj) <- "best cluster"
semisup15_markers <- FindAllMarkers(sobj )

# Filter markers for significant genes (e.g., adjusted p-value < 0.05 and log fold change > 0.25)
significant_markers <- semisup15_markers %>% 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

#VlnPlot(sobj, features = c("gene_of_interest"), group.by = "ident")

top10 <- significant_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sobj, features = top10$gene) + NoLegend()



Idents(sojb) <- "seurat_cluster.0.3"
DimPlot(sobj, reduction="umap", group.by='seurat_cluster.0.3')





############################
#


# 
# # cell counts by celltype as bar plot
# Idents(sobj2_tumor) <- "Sample.Label"
# #barplot(df_meta2, col=cols, las=2, cex.names=0.75, main="Cell counts by Sample")
# # by fov
# 
# df_meta <- sobj_tumor@meta.data
# df_meta2 <- sobj2_tumor@meta.data
# df_grp  <- df_meta %>%
#   group_by(Sample.Label) %>%
#   summarise(count = n())
# 
# 
# df_meta$count <- 1
# # Create the barplot
# ggplot(df_grp, aes(x = Sample.Label, y=count)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "Counts by Sample Grouped by Slide",
#        x = "Celltype",
#        y = "Count") +
#   theme_minimal()
# 
# table(df_meta$Sample.Label)
# 
# 
# 
# feature_names 
# # Filter feature names that begin with "Neg"
# neg_feature_names <- feature_names[grep("^Neg", feature_names)]
# print(neg_feature_names)
# 
# # Get row indices for neg_feature_names
# neg_indices <- which(feature_names %in% neg_feature_names)
# non_neg_indices <- which(!(feature_names %in% neg_feature_names))
# # Create sparse matrices for neg_feature_names and others
# neg_matrix <- counts_matrix[neg_indices, ]
# non_neg_matrix <- counts_matrix[non_neg_indices, ]
# 
# # Verify the results
# print(dim(neg_matrix))       # Dimensions of the matrix with neg_feature_names
# print(dim(non_neg_matrix))   # Dimensions of the matrix without neg_feature_names
# 
# cell_type_list <- sobj_CPA$Line
# group_averages <- calculate_group_averages(non_neg_matrix, cell_type_list)
# group_averages_matrix <- as.matrix(group_averages)
# 
# write.csv(as.data.frame(group_averages_matrix), 
#           "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_Tumors.csv",
#           row.names = TRUE)



#####################################
# PATH Annotations -> metadata
# needs to be done one slide at a time. 
####################################

# s_patient2 <- subset(sobj, subset = PatientID == 2 )
# df_meta <- s_patient2@meta.data
# plot_grouped_barplot(df_meta, cell.type='NS_Insitu_Celltype.level1', group.by='Sample.Label', title=" ")




sobj1 <- subset(sobj, subset = Run_Tissue_name == "PSR-01")
sobj2 <- subset(sobj, subset = Run_Tissue_name == "PSR-02")
sobj_org <- subset(sobj, subset = Run_Tissue_name == "PSR-Organoids")

# Dump cell coords
df_coords_psr1 <- sobj1@meta.data[,c("cell_id", "CenterX_global_px", "CenterY_global_px")]
head(df_coords_psr1)
df_coords_psr2 <- sobj2@meta.data[,c("cell_id", "CenterX_global_px", "CenterY_global_px")]
head(df_coords_psr2)
# write to csv
write.csv(df_coords_psr1, file = paste0(results_dir, "/df_coords_psr1.csv"))
write.csv(df_coords_psr2, file = paste0(results_dir, "/df_coords_psr2.csv"))

# get pixel class from mask files
path_mask_dir <- paste0("/Volumes/T7Shield/PSR-GEN-057/images/masks")
# no. pkg:tiff: Sorry, can not handle images with 64-bit samples.  Try saving these as lower bits. 
psr1_mask_img <- image_read(paste0(path_mask_dir, "/RNA_PSR01_pathology_mask_uint8.tif")) %>%
  image_convert(type = 'TrueColorAlpha') %>%
  #image_rotate(180) %>%
  image_flip() %>%
  # only black is transparent
  image_transparent("black", fuzz = 0) 
#image_modulate(brightness = 100, saturation = 100, hue = 100) 

raster_data <- as.raster(psr1_mask_img)
# Convert raster to dataframe for ggplot.  looks plausible but doesnt work
# raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)
# head(raster_df)
# colnames(raster_df) <- c("x", "y", "value")


# adding my own FOV isn't working; how to attach fov to sobj1? seems liks obj@images@new_fov <- fov
#sobj1 <- create_fov(sobj1, x_col="CenterX_global_px_new", y_col="CenterY_global_px_new", assay="SCT", name="PSRO1_small") 
#sobj1@images$image

#print(paste("dims of mask image: ", dim(psr1_mask_img))) # 1, 7K, 9K

df_annotated_xy <- get_pixel_class_df(psr1_mask_img)
df_annotated_xy <- as.data.frame(df_annotated_xy)
# Convert list of triples to a dataframe
head(df_annotated_xy)
nrow(df_annotated_xy)
p2 <- ggplot(df_annotated_xy, aes(x=x, y=y)) + geom_point()
p2

# ggplot both dataframes overlapping w cowplot
# cowplot supports overlapping images, and grid systems
#library(cowplot)

#p1 <- ggplot(df_coords_psr1, aes(x=CenterX_global_px_new, y=CenterY_global_px_new)) + geom_point()

################################
# see how close our rescale is, bringing points down to our image size
################################
image_width <- image_info(psr1_mask_img)$width
image_height <- image_info(psr1_mask_img)$height

# trial and error here to get the right offset.
offset_x <- 600 + 21 # adding to this moves stuff left
offset_y <- 200 - 21 # subtracting from this moves stuff up

putative_big_width <- image_width * 16  # 114784 # image_width * 16 # I got close with 114000 # x axis: larger number -> stretch more to right.
putative_big_height <- image_height * 16 # 144560 #  # I got close with 145000 

# new coords, rescaled to our tiff size (best we've got)
head(df_coords_psr1)
df_coords_psr1$CenterX_global_px_new <- round(df_coords_psr1$CenterX_global_px * ( image_width / putative_big_width)) + offset_x
df_coords_psr1$CenterY_global_px_new <- round(df_coords_psr1$CenterY_global_px * (image_height / putative_big_height )) + offset_y

sobj1 <- AddMetaData(sobj1, df_coords_psr1[,c("CenterX_global_px_new", "CenterY_global_px_new")])

df_annotated_xy <- df_annotated_xy %>% 
  mutate(path_annot_name = case_when(
    path_annot_id == 2 ~ "Adeno",
    path_annot_id == 3 ~ "Lymph-Inf",
    path_annot_id == 5 ~ "Tumor-SQ",
    TRUE ~ "Unknown"
  ))

p_path_annot_psr1 <- ggplot() +
  geom_point(data = df_coords_psr1, aes(x=CenterX_global_px_new, y=CenterY_global_px_new), color = "grey") +
  geom_point(data = df_annotated_xy, aes(x = x, y = y, color = path_annot_name)) + #color = c("blue", "lightblue", "pink")) +
  # this works on factors. 
  scale_color_manual(values = c("Adeno" = "lightblue", "Lymph-Inf" = "pink", "Tumor-SQ" = "blue")) + 
  labs(title = "Pathology Annotations over cells detected in PSR01",
       x = "", y = "")

p_path_annot_psr1

ggsave(paste0(results_dir, "/PSR01_whole_slide_path_annot_overlay.png"), p_path_annot_psr1, width = 10, height = 10, dpi = 300)


# FOUND scale of points is 16 times greater than the exported tiff image size. 
# is higher resolution an even mutliple of lower resolution?  
#114000 / image_width 
image_width * 16 # 16 found by trying different sizes.  offset needed, too.
image_height * 16

# TODO: move this to unittest
# Example data frame for testing 
# df <- data.frame(
#   x = c(1380, 1380, 1380),
#   y = c(8190, 8191, 8192),
#   path_annot_id = c(3, 3, 3)
# )
# 
# find_annotation (df_annotated_xy, 4, 4, tolerance=0)
# find_annotation (df_annotated_xy, 1379, 8194, tolerance=2)
# annotations <- find_annotation (df_annotated_xy, 1379, 8194, tolerance=1)
# print(annotations)

# must give in x and y, test with a couple points
#annotations <- find_annotation (df_annotated_xy, df_metadata$CenterX_global_px_new[2], df_metadata$CenterY_global_px_new[2])
#print(annotations)
#df_metadata$CenterX_global_px_new[2]
#df_metadata$CenterY_global_px_new[2]

# now sapply this to metadata on sboj
# this puts an empty list in each cell of the dataframe
df_metadata <- sobj1@meta.data

df_metadata$path_annot_id <- sapply(1:nrow(df_metadata), function(i) {
  find_annotation(as.data.frame(df_annotated_xy), df_metadata$CenterX_global_px_new[i], df_metadata$CenterY_global_px_new[i], tolerance=1)
})

head(df_metadata$path_annot_id)
table(df_metadata$path_annot_id)

df_metadata <- df_metadata %>%
  mutate( 
    path_annot_name = case_when(
      path_annot_id == 1 ~ "desmoplasia", # red 
      path_annot_id == 2 ~ "adeno", # lt blue
      path_annot_id == 3 ~ "lymph_inf", # pink
      path_annot_id == 5 ~ "tumor_sq", # blue 
      path_annot_id == 6 ~ "tumor_inf", # green
      path_annot_id == 9 ~ "SCLC", # yellow
      TRUE ~ "NA"
    )
  )

df_meta_to_add <- df_metadata %>%
  dplyr::select(path_annot_id, path_annot_name)

sobj1 <- AddMetaData(sobj1, df_meta_to_add)
head(sobj1@meta.data$path_annot_id)
saveRDS(sobj1, paste0(results_dir, "sobj1_path_annot.rds"))

sobj1 <- readRDS(paste0(results_dir, "sobj1_path_annot.rds"))

# Add path_annot_id and path_annot_name as metadata to sobj 
sobj <- AddMetaData(sobj, df_metadata[,c("path_annot_id", "path_annot_name")])

# Fill in all the NAS with 0
# dim(sobj)
# path_df <- sobj@meta.data[c("path_annot_id", "path_annot_name")]
# path_df[is.na(path_df)] <- 0
# table(path_df$path_annot_id)
# sobj <- AddMetaData(sobj, path_df)
# table(sobj@meta.data$path_annot_id) 

# TODO: Repeat for PSR02
psr2_mask_img <- image_read(paste0(path_mask_dir, "/RNA_PSR02_pathology_mask_uint8.tif"))



# save metatadata for use by giotto, also
df_metadata <- sobj@meta.data
write.csv(df_metadata, file = paste0(results_dir, "/full_metadata.csv"))

#############################################
# PSR 02

# Dump cell coords
# df_coords_psr1 <- sobj1@meta.data[,c("cell_id", "CenterX_global_px", "CenterY_global_px")]
# head(df_coords_psr1)
# df_coords_psr2 <- sobj2@meta.data[,c("cell_id", "CenterX_global_px", "CenterY_global_px")]
# head(df_coords_psr2)
# write to csv
#write.csv(df_coords_psr1, file = paste0(results_dir, "/df_coords_psr1.csv"))
#write.csv(df_coords_psr2, file = paste0(results_dir, "/df_coords_psr2.csv"))

# get pixel class from mask files
#path_mask_dir <- paste0(root_dir, '/images/masks')
# no. pkg:tiff: Sorry, can not handle images with 64-bit samples.  Try saving these as lower bits. 
psr2_mask_img <- image_read(paste0(path_mask_dir, "/RNA_PSR02_pathology_mask_uint8.tif")) %>%
  image_convert(type = 'TrueColorAlpha') %>%
  #image_rotate(180) %>%
  image_flip() %>%
  # only black is transparent
  image_transparent("black", fuzz = 0) 
#image_modulate(brightness = 100, saturation = 100, hue = 100) 

raster_data <- as.raster(psr2_mask_img)

df_annotated_xy2 <- get_pixel_class_df(psr2_mask_img)
df_annotated_xy2 <- as.data.frame(df_annotated_xy2)
# Convert list of triples to a dataframe
head(df_annotated_xy2)
nrow(df_annotated_xy2)
p2 <- ggplot(df_annotated_xy2, aes(x=x, y=y)) + geom_point()
p2

################################
# see how close our rescale is, bringing points down to our image size
################################
image_width <- image_info(psr2_mask_img)$width
image_height <- image_info(psr2_mask_img)$height

# trial and error here to get the right offset.
offset_x <- -280 # was -200 # adding to this moves stuff left
offset_y <- 440 # was 400 subtracting from this moves stuff up

putative_big_width <- image_width * 16  # 114784 # image_width * 16 # I got close with 114000 # x axis: larger number -> stretch more to right.
putative_big_height <- image_height * 16 # 144560 #  # I got close with 145000 

# new coords, rescaled to our tiff size (best we've got)
head(df_coords_psr2)
df_coords_psr2$CenterX_global_px_new <- round(df_coords_psr2$CenterX_global_px * ( image_width / putative_big_width)) + offset_x
df_coords_psr2$CenterY_global_px_new <- round(df_coords_psr2$CenterY_global_px * (image_height / putative_big_height )) + offset_y

sobj2 <- AddMetaData(sobj2, df_coords_psr2[,c("CenterX_global_px_new", "CenterY_global_px_new")])

df_annotated_xy2 <- df_annotated_xy2 %>%
  mutate(path_annot_name = case_when(
    path_annot_id == 1 ~ "desmoplasia", # red 
    path_annot_id == 2 ~ "adeno", # lt blue
    path_annot_id == 3 ~ "lymph_inf", # pink
    path_annot_id == 5 ~ "tumor_sq", # blue 
    path_annot_id == 6 ~ "tumor_inf", # green
    path_annot_id == 9 ~ "SCLC", # yellow
    TRUE ~ "Unknown"
  ))

annot_colors <- c("red","lightblue","pink", "blue", "green", "yellow")
#names(annot_colors) <- c("desmoplasia","adeno","lymph_inf", "tumor_sq", "tumor_inf", "SCLC")
annot_colors <- setNames(annot_colors, c("desmoplasia","adeno","lymph_inf", "tumor_sq", "tumor_inf", "SCLC"))

p_path_annot_psr2 <- ggplot() +
  geom_point(data = df_coords_psr2, aes(x=CenterX_global_px_new, y=CenterY_global_px_new), color = "grey") +
  geom_point(data = df_annotated_xy2, aes(x = x, y = y, color = path_annot_name)) + #color = c("blue", "lightblue", "pink")) +
  # this works on factors. 
  scale_color_manual(values = annot_colors) + 
  coord_fixed(ratio = 1) + 
  labs(title = "Pathology Annotations over cells detected in PSR02",
       x = "", y = "")

p_path_annot_psr2

ggsave(paste0(results_dir, results_sub_wholeassay, "/PSR02_whole_slide_path_annot_overlay.png"), p_path_annot_psr2, width = 10, height = 10, dpi = 300)

# now sapply this to metadata on sboj
# this puts an empty list in each cell of the dataframe
df_metadata <- sobj2@meta.data

df_metadata$path_annot_id <- sapply(1:nrow(df_metadata), function(i) {
 find_annotation(as.data.frame(df_annotated_xy2), df_metadata$CenterX_global_px_new[i], df_metadata$CenterY_global_px_new[i], tolerance=1)
})

# df1 <- df_metadata[,c("cell_id","CenterX_global_px_new", "CenterY_global_px_new")]
# head(df1)
# dim(df1)
# tolerance <- 1
# df_result <- df1 %>%
#   rowwise() %>%
#   mutate(
#     max_path_annot_id = max(
#       df_annotated_xy2$path_annot_id[abs(df1$CenterX_global_px_new - df_annotated_xy2$x) <= tolerance & 
#                           abs(df1$CenterY_global_px_new - df_annotated_xy2$y) <= tolerance], 
#       na.rm = FALSE  # In case no match is found
#     )
#   ) %>%
#   ungroup()

# started 6pm

head(df_metadata$path_annot_id)
table(df_metadata$path_annot_id)

df_metadata <- df_metadata %>%
  mutate( 
    path_annot_name = case_when(
      path_annot_id == 1 ~ "desmoplasia", # red 
      path_annot_id == 2 ~ "adeno", # lt blue
      path_annot_id == 3 ~ "lymph_inf", # pink
      path_annot_id == 5 ~ "tumor_sq", # blue 
      path_annot_id == 6 ~ "tumor_inf", # green
      path_annot_id == 9 ~ "SCLC", # yellow
      TRUE ~ "NA"
    )
  )

df_meta_to_add <- df_metadata %>%
  dplyr::select(path_annot_id, path_annot_name)

sobj2 <- AddMetaData(sobj2, df_meta_to_add)
head(sobj2@meta.data$path_annot_id)
saveRDS(sobj2, paste0(results_dir, "/sobj2_path_annot.rds"))

sobj2 <- readRDS(paste0(results_dir, "/sobj2_path_annot.rds"))
df_metadata <- sobj2@meta.data
table(df_metadata$path_annot_name)

table(sobj@meta.data$path_annot_name)

# Add path_annot_id and path_annot_name as metadata to sobj 
sobj <- AddMetaData(sobj, df_metadata[,c("path_annot_id", "path_annot_name")])

# Fill in all the NAS with 0
dim(sobj)
path_df <- sobj@meta.data[c("path_annot_id", "path_annot_name")]
path_df[is.na(path_df)] <- 0
table(path_df$path_annot_id)
sobj <- AddMetaData(sobj, path_df)
table(sobj1@meta.data$path_annot_name) 

# Repeat for PSR02
#psr2_mask_img <- image_read(paste0(path_mask_dir, "/RNA_PSR02_pathology_mask_uint8.tif"))

# save metatadata for use by giotto, also
df_metadata <- sobj@meta.data
write.csv(df_metadata, file = paste0(saved_rds_dir, "/full_metadata.csv"))



saveRDS(sobj, paste0(saved_rds_dir, "sobj_path_annot.RDS"))

# PSR01
offset_x <- 600 + 21 # adding to this moves stuff left
offset_y <- 200 - 21 # subtracting from this moves stuff up

df_meta_to_add <- add_whole_slide_annots_coords(sobj1, 
                                       paste0(path_mask_dir, "/RNA_PSR01_pathology_mask_uint8.tif"),
                                       offset_x = 621, 
                                       offset_y = 179, 
                                       "PSR01", 
                                       results_dir) 

# PSR02
offset_x <- -280 # was -200 # adding to this moves stuff left
offset_y <- 440

df_meta_to_add <- add_whole_slide_annots_coords(sobj2, 
                                       paste0(path_mask_dir, "/RNA_PSR02_pathology_mask_uint8.tif"),
                                       offset_x = -280, 
                                       offset_y = 440, 
                                       "PSR02", 
                                       results_dir) 

head(sobj2@meta.data$path_annot_id)


saveRDS(sobj, paste0(results_dir, "sobj_path_annot.RDS"))



saveRDS(sobj, file = paste0(saved_rds_dir, "/latest_sobj_30Sep24.rds"))
R.version.string


################################
# Pathway analysis by paired patient
################




################################
# plot mor gene differences
################  

# e.g. MET
patient_num = 3
DefaultAssay(sobj) <- "SCT"
s_patient3_pre <- subset(sobj_tumor, subset = PatientID == patient_num & TimePoint == "pre")
s_patient3_post <- subset(sobj_tumor, subset = PatientID == patient_num & TimePoint == "post")
dim(s_patient3_pre)
dim(s_patient3_post)

patient_num = 4
s_patient4 <- subset(sobj_tumor, subset = PatientID == patient_num)
dim(s_patient4)

patient_num = 6
s_patient6 <- subset(sobj_tumor, subset = PatientID == patient_num)
dim(s_patient6)

# Get the gene expression data for MET, but break out by TimePoint
gene_name = "MET"
gene_names = c("MET","EGFR", "KRT5")
gene_names = list(df_special_genes$Display_Name)
# gene_names = c("CCND1",  "CD9",    "EGFR",   "EPCAM",  "H2AZ1" , "KIT" ,   "KRT20",  "KRT5",   "KRT6",   "KRT7",   "MET",     
#                "RB1",    "S100B" ,
#                 "SRC",    "STAT3",  "TP53",   "ERBB2",  "ERBB3",  "NCAM1",  "FGFR1",  "MYC" ,   "HGF" ,   "JUP",    "BRAF" , 
#                "KRAS" ,  "TWIST1"
# )
# what is the default layer? data.
# gene_data_patient3_pre <- FetchData(s_patient3_pre, gene_name, layer = "data")
# gene_data_patient3_pre$TimePoint <- "pre"
# head(gene_data_patient3_pre)
# gene_data_patient3_post <- FetchData(s_patient3_post, gene_name, layer = "data")
# gene_data_patient3_post$TimePoint <- "post"
# 
# # metadata to add
# df_meta <- sobj@meta.data[c('TimePoint','Sample.ID','Sample.Label', 'PatientID')]
# head(df_meta)

# # combine with pre/post timepoint 
# gene_data_patient3 <- rbind(gene_data_patient3_pre, gene_data_patient3_post)
# # boxplot to compare TimePoint pre vs post
# head(gene_data_patient3)
# # join
# gene_data_patient3b <- merge(gene_data_patient3, df_meta, how="inner", by="row.names")
# head(gene_data_patient3b)

# test
#norm_data <- FetchData(sobj_CT_tumor_4, gene_names, layer = "data")
#df_meta <- sobj_CT_tumor_4@meta.data[c('TimePoint','Sample.ID','Sample.Label', 'PatientID')]
#gene_data <- merge(norm_data, df_meta, how="inner", by="row.names")


gene_data1 <- fetch_gene_data(sobj_PT_tumor4, gene_names)
head(gene_data1)
gene_data2 <- gene_data1 %>%
  dplyr::filter(PatientID %in% c(3,4,6)) 


# Assuming your dataset is called df
ggplot(gene_data2, aes(x = TimePoint, y = MET)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "MET values by TimePoint", x = "TimePoint", y = "MET") 

# Subset the data by TimePoint
pre <- gene_data_patient3$MET[gene_data_patient3$TimePoint == "pre"]
post <- gene_data_patient3$MET[gene_data_patient3$TimePoint == "post"]
# Perform paired t-test
t_test_result <- t.test(pre, post, paired = FALSE)
# View the p-value from the t-test
t_test_result$p.value

table(gene_data2$Sample.Label)
# Create a boxplot with 6 boxes based on the 'Group' variable
Group= gene_data2$Sample.Label
# control display order
gene_data2$Group <- factor(gene_data2$Sample.Label, levels = c("3-lymph node-pre ( 3_5 )", 
                                                        "3-liver-post ( 3_6 )", 
                                                        "4-liver-pre ( 4_6 )", 
                                                        "4-liver-post ( 4_1 )",
                                                        "6-brain-pre ( 6_3 )",
                                                        "6-lymph node-post ( 6_5 )"))
gene_data2$Group <- factor(gene_data2$TimePoint, levels = c("pre", "post"))


marker = "MET"
# change next line too. 
ggplot(gene_data2, aes(x = Group, y = MET, group=Sample.Label, fill=Group)) +
  geom_boxplot() +
  #theme_minimal() +
  #scale_fill_manual(values = c("Group1" = "#FF0000", "Group2" = "#00FF00", "Group3" = "#0000FF", "Group4" = "#FFFF00")) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +
  labs(title = paste(marker,"values by Group"), x = "Group", y = marker) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(results_dir, "/boxplot_Pt4_", marker, ".png"), last_plot(), width = 4, height = 4, dpi = 300)
# loop over


# how to plot stack of multiple genes
# Create the faceted boxplots
#library(tidyr)

file_suffix <- "panckTumor"
file_suffix <- "celltypeTumor"

# Reshape the dataframe from wide to long format
gene_data_long <- gene_data2 %>%
  pivot_longer(cols = c(MET, EGFR, KRT5), names_to = "Gene", values_to = "Expression")

# Now you can plot using facet_wrap
ggplot(gene_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  labs(title = "Boxplots for Multiple Genes", x = "Group", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# adeno
gene_names <- c("EGFR", "SLPI","KRT5","EPCAM","KRT7","KRT20","MET")
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("NSCLC markers in Tumor Cells ", file_suffix))
p 
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_NSCLC_", gene_name, ".png"), p, width = 10, height = 40, dpi = 300)

# Oncogenes
gene_names <- c("EGFR", "ERBB2", "JUN","JUNB", "KRAS", "RELA", "SRC","YES1")
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle(paste0("Oncogene markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Oncogenes.png"), p, dpi = 300)


# Small Cell # too low
gene_names <- c("YES1")
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle("Small Cell markers in Tumor Cells")
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_SmallCell.png"), p, dpi = 300)

# Squamous # too low expression to tell 
gene_names <- c("KRT5", "KRT6A/B/C","KRT13","KRT14") # TP63
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names)+
  ggtitle(paste0("Squamous Cell markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Squamous.png"), p, dpi = 300)

# NE, TP53+ RB1- c be NE. Too low 
# NCAM1 in panel list but not found
gene_names <- c("TP53", "RB1","S100B") # , "NCAM1") # , "NCAM", "SYP" )
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("NE markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_NE1.png"), p,  dpi = 300)

# e.g. super low
gene_data1$RB1


# Plasticity
gene_names <- c("EGFR", "SLPI","KRT5","EPCAM","KRT7","KRT20", "TP53","RB1","MET")
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("Plasticity markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Plasticiy_", gene_name, ".png"), p, width = 10, height = 40, dpi = 300)


# Gentles adeno
# NE, TP53+ RB1- c be NE. Too low 
gene_names <- gentles_adeno
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("Gentles Adeno markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Gentles_Adeno.png"), p,  dpi = 300)
# Gentles sq
gene_names <- gentles_sq
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("Gentles Squamous markers in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_Gentles_Sq.png"), p,  dpi = 300)


# have boxes
gene_names <- c("SLPI","EPCAM","KRT7","MET", "JUNB")
gene_data1 <- fetch_gene_data(sobj_panck_tumor, gene_names)
p <- plot_multiple_gene_boxplots(gene_data1, gene_names) +
  ggtitle(paste0("Highest Signal LC markers of interest in Tumor Cells ", file_suffix))
p
ggsave(paste0(results_dir, "/boxplot_",file_suffix,"_MostSignal.png"), p,  dpi = 300)



# Violin Plot of multiple genes
# Using an external vector in selections was deprecated in tidyselect 1.1.0.
# ℹ Please use `all_of()` or `any_of()` instead.
# # Was:
# data %>% select(gene_list)
# 
# # Now:
# data %>% select(all_of(gene_list))
plot_multiple_gene_boxplots <- function(df_gene_data, gene_list, group.by="Sample.Label"){
  
  if (group.by == "Sample.Label") {
  df_gene_data$Group <- factor(df_gene_data$Sample.Label, 
                               levels = c("3-lymph node-pre ( 3_5 )", 
                                          "3-liver-post ( 3_6 )", 
                                          "4-liver-pre ( 4_6 )", 
                                          "4-liver-post ( 4_1 )",
                                          "6-brain-pre ( 6_3 )",
                                          "6-lymph node-post ( 6_5 )"))
  }
  else {
    df_gene_data$Group <- factor(df_gene_data[group.by])
  }
  
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


gene_data_patient4 <- FetchData(s_patient4, gene_name)
gene_data_patient6 <- FetchData(s_patient6, gene_name)
gene_data_patient3


# Also do this with tumor cell types. 


# weird stuff 
# dim(df_gene_data)
# df_gene2 <- df_gene_data %>% as.data.frame()
# type(df_gene2)
# df_gene_data1 <- df_gene_data$Display_Name
# df_gene_data2 <- df_gene_data$Gene_Name_s_
# df_gene_data_new <- data.frame(df_gene_data1, df_gene_data2)
# colnames(df_gene_data_new) <- c("Display_Name", "descr")
# 
genes_annot <- df_pre %>%
  dplyr::left_join(df_gene_data_new, by = c("gene" = "Display_Name"))

head(genes_annot)
# Get top 30 abs value of log2FC
top_genes <- genes_annot %>%
  top_n(30, wt = abs(avg_log2FC))

# Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 4), lab_fill = "zscore")

#WaterfallPlot(genes_annot, f = sobj_path_tumor$Timepoint, ident.1 = "pre", ident.2 = "post", top.n = 20)
# Create a waterfall plot with ggplot2
ggplot(top_genes, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_bar(stat = "identity") +   # Create the bar plot
  coord_flip() +                  # Flip the coordinates to make it horizontal
  geom_text(aes(label = descr),    # Add descriptions as labels
                  size = 3, hjust = 0, vjust = 0.5) + 
  scale_fill_manual(values = c("red", "blue"), guide = "none") +  # Red for negative, blue for positive
  theme_minimal() + 
  labs(title = "Waterfall Plot of Top DE Genes in Pre-Treatment Tumor Samples",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for clarity

write.csv(genes_annot, paste0(results_dir, "/top_genes_pre_pt2.csv"))

##########################
# Gene Set Analysis Setup
#########################
# parent examples: lookup up https://amigo.geneontology.org/amigo/dd_browse
# immune_system_process
# cellular process GO:0050875, cell cycle process, cell cycle, regulation of cellular process
# ALK: check neurogenesis, cell proliferation, mapk signaling pathway, protein tyrosine kinase activity,
#     ATP binding GO:0140657, plasma membrane (cellular component). 
# 
options(max.print = 12, spe = "human")
getOption("spe")
#getGODatabase(parent = "cellular_process", root = "GO:0008150", spe = "human")
# biological_process = GO:0008150 worked, but takes forever. 
# Try these subsets for cancer
# "Regulation of immune response" (GO:0050776) - 
# "Negative regulation of apoptotic process" (GO:0043066): 
# "Cellular response to stress" (GO:0033554)
# "Regulation of cell proliferation" (GO:0042127):
# "Angiogenesis" (GO:0001525):
# "Extracellular matrix disassembly" (GO:0022617): n/a

# which do I have so far? 
l1 <- list_go_matrices(sobj)
l1
sobj <- add_go_matrix(sobj, go_parent="immune_system_process") # 118
sobj <- add_go_matrix(sobj, go_parent="GO:0050776") # 106
sobj <- add_go_matrix(sobj, go_parent="GO:0043066") # 18
sobj <- add_go_matrix(sobj, go_parent="GO:0033554") # 11 
sobj <- add_go_matrix(sobj, go_parent="GO:0042127") # 87
sobj <- add_go_matrix(sobj, go_parent="GO:0001525") # 3

# sobj with GSEA matrices
saveRDS(sobj, paste0(saved_rds_dir, "/latest_sobj_3Oct24.rds"))

# lets try just getting the matrix
dim(sobj@misc$AUCell$GO[['GO:0050776']])

matr <- sobj@misc$AUCell$GO[['GO:0050776']]
dim(matr)
matr <- RenameGO(matr)

WaterfallPlot(matr, f = sobj$CellType_Labeled.2, 
              ident.1 = "Fibroblast", ident.2 = "Macrophage",
              len.threshold = 1)

# re-do subsets to get this GSEA
sobj_tumor <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") &
                                 CellType_Labeled.1 == "Tumor" &
                                 TimePoint %in% c("pre", "post"))
sobj_tumor_pairs <- subset(sobj_tumor, subset = (PatientID %in% c(3,4,6) & Sample.ID != "ORG_4_1")
                           | Run_Tissue_name == "PSR-CPA")

sobj_pt3_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 3)
sobj_pt4_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 4)
sobj_pt6_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 6)

Idents(sobj_pt3_tumor) <- "CellType_Labeled"
table(sobj_pt3_tumor$CellType_Labeled)
Idents(sobj_pt4_tumor) <- "CellType_Labeled"
table(sobj_pt4_tumor$CellType_Labeled)
Idents(sobj_pt6_tumor) <- "CellType_Labeled"
table(sobj_pt6_tumor$CellType_Labeled)


run_gsea <- function(sobj, group.by="Sample.Label2", ident1, ident2, go_id, go_name, th=2, patientid, results_dir){
  matr <- sobj@misc$AUCell$GO[[go_id]]
  matr <- RenameGO(matr)
  print(paste0("matr has ", dim(matr)[1], " rows."))
  Idents(sobj) <- group.by
  p <- WaterfallPlot(matr, f = sobj$Sample.Label2,  # change to group.by?
              ident.1 = ident1, ident.2 = ident2,
              len.threshold = th) +
    theme(text = element_text(size = rel(2))) +
  labs(title = paste0("GSEA Plot of ", go_name, " in Pre-Post-Treatment Tumor Samples, Pt", patientid)) 
  
  p
  ggsave(paste0(results_dir, "/GSEA_pre_post_tumor_pt",patientid, "_", sanitize_name(go_name), ".png"), 
         p, height=10, width=20, dpi=300)
  return (p)
}

p <- run_gsea(sobj_pt3_tumor, group.by="Sample.Label2", ident1 = "CUTO-59-lymph node-pre (3_5)", ident2 = "CUTO-59-liver-post (3_6)",
         go_id = "GO:0050776", go_name = "Regulation of Immune Response", th=2, patientid=3, results_dir=results_dir)
p
p <- run_gsea(sobj_pt4_tumor, group.by="Sample.Label2", ident1 = "220904-liver-pre (4_6)", ident2 = "220904-liver-post (4_1)",
         go_id = "GO:0050776", go_name = "Regulation of Immune Response", th=2, patientid=4, results_dir=results_dir)
p
p <- run_gsea(sobj_pt6_tumor, group.by="Sample.Label2", ident1 = "220901-brain-pre (6_3)", ident2 = "220901-lymph node-post (6_5)",
         go_id = "GO:0050776", go_name = "Regulation of Immune Response", th=4, patientid=6, results_dir=results_dir)
p

############# Next GO term "Negative regulation of apoptotic process" (GO:0043066): 
go_id <- "GO:0043066"
dim(sobj_pt3_tumor@misc$AUCell$GO[[go_id]])
p <- run_gsea(sobj_pt3_tumor, group.by="Sample.Label2", ident1 = "CUTO-59-lymph node-pre (3_5)", ident2 = "CUTO-59-liver-post (3_6)",
              go_id = "GO:0043066", go_name = "Negative reg of apoptotic process", th=2, patientid=3, results_dir=results_dir)
p
p <- run_gsea(sobj_pt4_tumor, group.by="Sample.Label2", ident1 = "220904-liver-pre (4_6)", ident2 = "220904-liver-post (4_1)",
              go_id = "GO:0043066", go_name = "Negative reg of apoptotic process", th=2, patientid=4, results_dir=results_dir)
p
p <- run_gsea(sobj_pt6_tumor, group.by="Sample.Label2", ident1 = "220901-brain-pre (6_3)", ident2 = "220901-lymph node-post (6_5)",
              go_id = "GO:0043066", go_name = "Negative reg of apoptotic process", th=3, patientid=6, results_dir=results_dir)
p

# "Cellular response to stress" (GO:0033554)
p <- run_gsea(sobj_pt3_tumor, group.by="Sample.Label2", ident1 = "CUTO-59-lymph node-pre (3_5)", ident2 = "CUTO-59-liver-post (3_6)",
              go_id = "GO:0033554", go_name = "Cellular response to stress", th=1, patientid=3, results_dir=results_dir)
p
p <- run_gsea(sobj_pt4_tumor, group.by="Sample.Label2", ident1 = "220904-liver-pre (4_6)", ident2 = "220904-liver-post (4_1)",
              go_id = "GO:0033554", go_name = "Cellular response to stress", th=1, patientid=4, results_dir=results_dir)
p
p <- run_gsea(sobj_pt6_tumor, group.by="Sample.Label2", ident1 = "220901-brain-pre (6_3)", ident2 = "220901-lymph node-post (6_5)",
              go_id = "GO:0033554", go_name = "Cellular response to stress", th=1, patientid=6, results_dir=results_dir)
p

#"Regulation of cell proliferation" (GO:0042127):
p <- run_gsea(sobj_pt3_tumor, group.by="Sample.Label2", ident1 = "CUTO-59-lymph node-pre (3_5)", ident2 = "CUTO-59-liver-post (3_6)",
              go_id = "GO:0042127", go_name = "Reg of cell proliferation", th=2, patientid=3, results_dir=results_dir)
p
p <- run_gsea(sobj_pt4_tumor, group.by="Sample.Label2", ident1 = "220904-liver-pre (4_6)", ident2 = "220904-liver-post (4_1)",
              go_id = "GO:0042127", go_name = "Reg of cell proliferation", th=2, patientid=4, results_dir=results_dir)
p
p <- run_gsea(sobj_pt6_tumor, group.by="Sample.Label2", ident1 = "220901-brain-pre (6_3)", ident2 = "220901-lymph node-post (6_5)",
              go_id = "GO:0042127", go_name = "Reg of cell proliferation", th=2, patientid=6, results_dir=results_dir)
p


# To run Waterfall gene FindMarkers


########### Pt 3
sobj_pt3_tumor <- GeneSetAnalysis(sobj_pt3_tumor, genesets = hall50$human)
matr <- sobj_pt3_tumor@misc$AUCell$genesets
WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "Mono CD14", ident.2 = "CD8 T cell")

table(sobj_pt3_tumor$Sample.Label2)
table(sobj_pt3_tumor@meta.data$CellType_Labeled.2)
genes <- VariableFeatures(sobj_pt3_tumor)[1:20]
WaterfallPlot(
  sobj_pt3_tumor, group.by = "Sample.Label2", features = genes,
  ident.1 = "CUTO-59-lymph node-pre (3_5)", ident.2 = "CUTO-59-liver-post (3_6)", length = "logFC",
  )
ggsave(paste0(results_dir, "/pt3_tumor_waterfall.png"), last_plot(), width = 30, height = 6, dpi = 300)

########### Pt 4
table(sobj_pt4_tumor$Sample.Label2)
#table(sobj_pt4_tumor@meta.data$CellType_Labeled.2)
genes <- VariableFeatures(sobj_pt4_tumor)[1:30]
WaterfallPlot(
  sobj_pt4_tumor, group.by = "Sample.Label2", features = genes,
  ident.1 = "220904-liver-pre (4_6)", ident.2 = "220904-liver-post (4_1)", length = "logFC")
ggsave(paste0(results_dir, "/pt4_tumor_waterfall.png"), last_plot(), width = 10, height = 6, dpi = 300)


table(sobj_pt6_tumor$Sample.Label2)
#table(sobj_pt4_tumor@meta.data$CellType_Labeled.2)
genes <- VariableFeatures(sobj_pt6_tumor)[1:20]
WaterfallPlot(
  sobj_pt6_tumor, group.by = "Sample.Label2", features = genes,
  ident.1 = "220901-brain-pre (6_3)", ident.2 = "220901-lymph node-post (6_5)", length = "logFC")
ggsave(paste0(results_dir, "/pt6_tumor_waterfall.png"), last_plot(), width = 20, height = 6, dpi = 300)



VolcanoPlot(sobj_path_tumor, group.by = "TimePoint", ident.1 = "pre", ident.2 = "post", min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay="SCT")



DoHeatmap(sobj_path_tumor, features = VariableFeatures(sobj_path_tumor)[1:50])
ggsave(paste0(results_dir, "/heatmap_Pt2_top50_heatmap.png"), last_plot())
FeaturePlot(s_new, features = gentles_adeno, ncol = 2)
RidgePlot(s_new, features = gentles_adeno, ncol = 2)

# volcano plot
# 2_3 vs 2_4

VolcanoPlot(s_new, group.by = "Sample.Label", ident.1 = "2_3", ident.2 = "2_4", min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay="SCT")





###################################################
# Crop and plot Celltypes
# CROP doesn't work because of different orientations of polygons and points
# Any molecules plot zooms out to full slide even if only one FOV is subset. 
##################################################
# subset to fov 24, add bg image, and ImageDimPlot + SpatialPlot


sobj1 <- subset(sobj, subset = Run_Tissue_name == "PSR-01")
df_metadata <- sobj1@meta.data
write.csv(df_metadata, file = paste0(results_dir, "/full_metadata1.csv"))


s_lung_10 <- subset(sobj, subset = fov == 10) # liver, contains some immune cells
s_lung_34 <- subset(sobj, subset = fov == 34) # LN

sobj1@images$PSR01$segmentation@bbox # max dimensions, origin top right. 
s_lung_34@images$PSR01$segmentation@bbox
s_lung_11@images$PSR01$segmentation@bbox
s_lung_24@images$PSR01$segmentation@bbox
s_lung_27@images$PSR01$segmentation@bbox
s_lung_26@images$PSR01$centroids@bbox # aligns with CenterX_global_px, CenterY_global_px coords, origin bottom left.


# check coords and bbox of whole slide
sobj1@images$PSR01$segmentation@bbox
ImageDimPlot(sobj34, fov="PSR01", )

#################
# LUNG + matching Organoid Cell Typing by Clustering
#################

s_lung_tumor <- subset(s_lung, subset = path_annot_id > 0)



s_lung <- subset(sobj1, subset = tissue == "lung" & fov %in% c(24, 26, 27))

# get coords of fov 24
df_meta <- sobj_spleen1@meta.data
df_meta <- df_meta %>%
  filter(Sample.ID == 31)
colnames(df_meta)
padding <- 5
x_min <- round(min(df_meta$CenterX_global_px)) - padding
y_min <- round(min(df_meta$CenterY_global_px)) - padding
x_max <- round(max(df_meta$CenterX_global_px)) + padding
y_max <- round(max(df_meta$CenterY_global_px)) + padding
# x_min <- round(min(df_meta$CenterX_local_px)) 
# y_min <- round(min(df_meta$CenterY_local_px))
# x_max <- round(max(df_meta$CenterX_local_px)) 
# y_max <- round(max(df_meta$CenterY_local_px)) 
print(paste(y_min, y_max, x_min, x_max)) #  "64632 68855 73854 78077"
# vs check FOV that was automatically resized with subsetting
head(s_lung@images$PSR01$centroids@coords) # like x: 74636 y: 68256 in range x: 3460-78072, y: 9817-68850
print(paste("x centroids range:",round(min(s_lung@images$PSR01$centroids@coords[,1])), max(s_lung@images$PSR01$centroids@coords[,1]))) # 3460 78072
print(paste("y centroids range:",round(min(s_lung@images$PSR01$centroids@coords[,2])), max(s_lung@images$PSR01$centroids@coords[,2]))) # 9817 68850
s_lung@images$PSR01$segmentation@bbox 
# min    max
# x 23472  98175
# y 70956 129975
# which is the right size for a FOV but with different origin? 

# create a Crop
# Note: cropped.coords should end up with a list of 3 named objects: centroids, segmentation, and molecules
# try by centroids or segmentation   
ImageDimPlot(s_lung_27, fov="PSR01", flip_xy=FALSE, axes=TRUE)
bbox <- s_lung_27@images$PSR01$centroids@bbox
bbox
# CROP
# coords: plot will crop from coords shown when plotting, "tissue" coords from GetTissueCoordinates. Are those the same?
sobj_spleen1@images[1]

cropped.coords <- Crop(s_lung_27[["PSR01"]], 
                       x = c(bbox['x','min'], bbox['x','max']), y = c(bbox['y','min'], bbox['y','max']), coords = "tissue")
# s/b same as this
cropped.coords <- Crop(sobj_spleen1[["Sabaawynewcore091320245"]], 
                       x = c(15108, 19160), y = c(9817,13919), coords = "tissue")

#my_fov <- CreateFOV(coords = coords_df, type = "centroids", assay="RNA", key="PSR-02")
sobj_spleen1@images$test <- cropped.coords
ImageDimPlot(sobj_spleen1, fov="test") 
ImageDimPlot(sobj_spleen1, fov="Sabaawynewcore091320245")

# Can I add the polygons?
##########################
# POLYGONS?
sobj_spleen1@images$Sabaawynewcore091320245$segmentation@bbox
sseg <- sobj_spleen1@images$Sabaawynewcore091320245$segmentation
class(sseg)




# earlier attempts will first object to the region not containing any centroids.  assume this *does* include centroids, but...
# The selected region does not contain any cell segmentations  
bbox_seg <- sobj_spleen1@images$Sabaawynewcore091320245$segmentation@bbox
bbox_seg
x_min <- bbox_seg[1,1]
x_max <- bbox_seg[1,2]
y_min <- bbox_seg[2,1]
y_max <- bbox_seg[2,2]

data <- c(x_min, y_min, x_max, y_max)



# get sample bbox
get_bbox_of_sample <- function(df_meta, fov_name, sample_id){

  df_meta <- df_meta %>%
    filter(Sample.ID == sample_id & Run_Tissue_name == fov_name)
  colnames(df_meta)
  padding <- 0
  x_min <- round(min(df_meta$CenterX_global_px)) - padding
  y_min <- round(min(df_meta$CenterY_global_px)) - padding
  x_max <- round(max(df_meta$CenterX_global_px)) + padding
  y_max <- round(max(df_meta$CenterY_global_px)) + padding

  data <- c(x_min, y_min, x_max, y_max)
  # Create the matrix and specify row and column names
  bbox_mtx <- matrix(data, nrow = 2, byrow = FALSE,
                      dimnames = list(c("x", "y"), c("min", "max")))
  print(bbox_mtx)   
 
  return (bbox_mtx)   
}

df_meta <- sobj@meta.data

bbox <- get_bbox_of_sample(df_meta, "Sabaawy new core 09/13/2024 5", 34)
ImageDimPlot(sobj, fov="Sabaawynewcore091320246", flip_xy=FALSE, axes=TRUE)
df_meta <- sobj@meta.data
df_meta <- df_meta %>%
  filter(fov == 1)
head(df_meta["c_1_1_1", c("CenterX_global_px", "CenterY_global_px", "Sample.ID")])

## bbox for sample 7



# no
cropped.coords <- Crop(sobj_spleen1[["Sabaawynewcore091320245"]], 
                       x = NULL, y = c(bbox_seg['y','min'], bbox_seg['y','max']), 
                       coords = c("plot", "tissue"))

# Also: Crop Error: Cannot remove default boundary 

# This returns labeld cells with x and y
GetTissueCoordinates(s_lung_27)
GetTissueCoordinates(sobj_spleen1)

# set a new field of view (fov)
s_lung_24[["fov24"]] <- cropped.coords

# visualize FOV using default settings (no cell boundaries)
p1 <- ImageDimPlot(vizgen.obj, fov = "hippo", axes = TRUE, size = 0.7, border.color = "white", cols = "polychrome",
                   coord.fixed = FALSE)


# Dim redux only lung
s_lung <- SCTransform(s_lung, assay = "Nanostring", verbose = TRUE)
s_lung_tumor <- SCTransform(s_lung_tumor, assay = "Nanostring", verbose = TRUE)

VariableFeatures(s_lung) <- VariableFeatures(s_lung)
VariableFeatures(s_lung_tumor) <- VariableFeatures(s_lung_tumor)

p1a <- VariableFeaturePlot(s_lung) + ggtitle("All lung cells")
top10 <- head(VariableFeatures(s_lung))
p2a <- LabelPoints(plot = p1a, points=top10, repel=TRUE)
#p1a + p2a

p1b <- VariableFeaturePlot(s_lung_tumor) + ggtitle("Annotated Tumor only")
top10 <- head(VariableFeatures(s_lung_tumor))
p2b <- LabelPoints(plot = p1b, points=top10, repel=TRUE) 
p2a + p2b
ggsave(paste0(results_dir, "VariableFeaturePlot_lung_all_v_tumor.png"), last_plot())

s_lung <- RunPCA(s_lung, verbose = FALSE, npcs = 20)
s_lung_tumor <- RunPCA(s_lung_tumor, verbose = FALSE, npcs = 20)
# What features account for the most variation within lung samples?
p_lung_dim_loadings <- VizDimLoadings(s_lung, dims = 1:2, reduction = "pca") + ggtitle("Largest sources of variation")
ggsave(paste0(results_dir, "VizDimLoadings_lung.png"), p_lung_dim_loadings)

p_lung_tumor_dim_loadings <- VizDimLoadings(s_lung_tumor, dims = 1:2, reduction = "pca") + ggtitle("Largest sources of variation tumor only")
p_lung_tumor_dim_loadings
ggsave(paste0(results_dir, "VizDimLoadings_lung.png"), p_lung_tumor_dim_loadings)


# Assess how PCA did with different npcs, 20, 30, 50 e.g. 
# ElbowPlot to decide best cutoff.  DimPlot colored by possible batch effect variables. e.g. slides should look integrated.
ElbowPlot(s_lung) # 14 PCS is a good cutoff or even 7,
ElbowPlot(s_lung_tumor) # 9 PCS is a good cutoff or even 5,
#s_lung <- RunPCA(s_lung, verbose = FALSE, npcs = 15)
s_lung <- RunPCA(s_lung, verbose = FALSE, npcs = 7, features=VariableFeatures(s_lung))
DimHeatmap(s_lung, dims = 1:7, cells = 500, balanced = TRUE)
s_lung_tumor <- RunPCA(s_lung_tumor, verbose = FALSE, npcs = 5, features=VariableFeatures(s_lung_tumor))
DimHeatmap(s_lung_tumor, dims = 1:5, cells = 500, balanced = TRUE)

# Check for batch effects?
DimPlot(s_lung, reduction = "pca", group.by = "Sample.Label", alpha=0.5, cols='glasbey')
DimPlot(s_lung, reduction = "pca", group.by = "PanCK.PT", alpha=0.5, cols='glasbey')
dim(s_lung@meta.data) # only 2400 cells doesn't support many dimensions
DimPlot(s_lung_tumor, reduction = "pca", group.by = "Sample.Label", alpha=0.5, cols='glasbey')
DimPlot(s_lung_tumor, reduction = "pca", group.by = "fov", alpha=0.5, cols='glasbey')
DimPlot(s_lung_tumor, reduction = "pca", group.by = "PanCK.PT", alpha=0.5, cols='glasbey')


# What features explain most of the lung variance?

#
s_lung <- FindNeighbors(s_lung, dims = 1:7)
# Computes nearest neighbor graph, SNN in sobj@graphs$SCT_nn, SCT_nn. This is used for clustering.
s_lung <- FindClusters(s_lung, resolution=0.3, cluster.name='seurat_cluster.0.3', verbose = TRUE) # try different values of resolution for more/less clust.
# lowest resolution (0.1 finds 4 communities, 0.3 finds 9 )
s_lung <- RunUMAP(s_lung, dims = 1:7)

# seruat_cluster.0.2 is now at attribute in metadata
p_lung_clusters <- DimPlot(s_lung, group.by='seurat_cluster.0.3')
p_lung_clusters
DimPlot(s_lung, group.by='fov')


s_lung_tumor <- FindNeighbors(s_lung_tumor, dims = 1:5)
s_lung_tumor <- FindClusters(s_lung_tumor, resolution=0.2,  cluster.name='seurat_cluster.0.2', verbose = TRUE)
s_lung_tumor <- RunUMAP(s_lung_tumor, dims = 1:5)
p_lung_tumor_clusters <- DimPlot(s_lung_tumor, group.by='seurat_cluster.0.2') + ggtitle("Lung Tumor only - cell clusters")

# group.by is same as setting Ident
DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5, group.by='seurat_cluster.0.3')
ggsave(paste0(results_dir, "DimPlot_lung_8Clusters.png"), last_plot())

DimPlot(s_lung_tumor, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(paste0(results_dir, "DimPlot_lung_tumorOnly_5Clusters.png"), last_plot())

# # Fine tuning UMAP, range on n_neighbors, min_dist. For lung, no real difference. 
# n_neighbors = c(10, 20, 30, 40)
# min_dist = c(0.1, 0.2, 0.3, 0.4)
# # These all look identical.
# s_lung <- RunUMAP(s_lung, dims = 1:10, min_dist = 0.1, n_neighbors = 30)
# DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5)
# s_lung <- RunUMAP(s_lung, dims = 1:10, min_dist = 0.2, n_neighbors = 30)
# DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5)
# s_lung <- RunUMAP(s_lung, dims = 1:10, min_dist = 0.3, n_neighbors = 30)
# DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5)
# s_lung <- RunUMAP(s_lung, dims = 1:10, min_dist = 1, n_neighbors = 30)
# DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5)


### Cell Typing 
Idents(s_lung) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" 
DimPlot(s_lung, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(paste0(results_dir, "DimPlot_lung_9Clusters_CellType.png"), last_plot())

Idents(s_lung_tumor) <- 'NS_Insitu_Celltype.level1'
DimPlot(s_lung_tumor, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("Lung Tumor Only NS Insitu Cell Types L1")
ggsave(paste0(results_dir, "DimPlot_lung_tumorOnly_Clusters_CellType.png"), last_plot())

Idents(s_lung) <- "seurat_cluster.0.6"
Idents(s_lung) <- "path_annot_name"
p_immune_molecules <- ImageDimPlot(s_lung, 
             fov="PSR01",
             cells = row.names(s_lung@meta.data)[which(s_lung@meta.data$fov == 24)],
             alpha = 0.5,
             molecules = c("PTPRC", "CD8A", "CD4", "CD68", "CD19", "FOXP3", "Mean.CD45"), #
             #"EPCAM", "KRT5", "KRT7", "Mean.PanCK")
             mols.size = 0.02,
             axes = FALSE) + 
  ggtitle("Immune Molecules in Lung - Bug: showing full slide, not just fov")
p_immune_molecules

p_tumor_molecules <- ImageDimPlot(s_lung, 
                                   fov="PSR01",
                                   cells = row.names(s_lung@meta.data)[which(s_lung@meta.data$fov == 24)],
                                   alpha = 0.5,
                                   molecules = c("EPCAM", "KRT5", "KRT7", "PanCK.PT"), #
                                   #"EPCAM", "KRT5", "KRT7", "Mean.PanCK")
                                   mols.size = 0.2,
                                   axes = FALSE) + 
  ggtitle("Tumor Molecules in Lung - Bug: showing full slide, not just fov")
p_tumor_molecules

# let's instead label the clusters ourselves
#s_lung <- RenameIdents(s_lung, c("0" = "Tumor", "1" = "Immune", "2" = "Immune", "3" = "Immune", "4" = "Immune", "5" = "Immune", "6" = "Immune", "7" = "Immune", "8" = "Immune"))

Idents(s_lung) <- "seurat_cluster.0.3"
FeaturePlot(s_lung, features = c("PTPRC", "CD3E", "CD8A", "CD4", "CD68", "CD19", "FOXP3", "Mean.CD45"))
FeaturePlot(s_lung, features = c("EPCAM", "KRT5", "KRT7", "Mean.PanCK"))
DimPlot(s_lung, group.by = "NS_Insitu_Celltype.level1")
DimPlot(s_lung, group.by = "Sample.Label")
DimPlot(s_lung, group.by = "fov")

ImageDimPlot()

ImageFeaturePlot(s_lung_24, fov="PSR01", dark.background=FALSE, size=5, features= "EPCAM")
                 , "Mean.CD3E", "Mean.CD8A", "Mean.CD4", "Mean.CD19", "Mean.FOXP3", "Mean.KRT5", "Mean.KRT7", "Mean.EPCAM"))

s_lung_23 <- subset(s_lung, subset = fov == "23")
ImageFeaturePlot(s_lung_23, fov="PSR01", dark.background=FALSE, size=3, features= c("CD45", "CD8A", "CD4"))

# Idents choices
# "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" 
# "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"
# nn_59b0955e.a82f.4075.acc4.dac765dd45d6_1_cluster_cluster_2c2d99e8.e433.468c.8f08.066b386a1742_1"  
# "RNA_nbclust_c5e1d3d8.0a17.443a.994b.9a20f297fd2c_1_clusters"  
# NS_Insitu_Celltype.level1
# how to set ident to the clusters from umap

DefaultBoundary(s_lung@images$PSR01) <- "centroids"
DefaultBoundary(s_lung@images$PSR01) <- "segmentation"
Idents(s_lung_tumor) <- "NS_Insitu_Celltype.level1"
s_lung_24 <- subset(s_lung, subset = fov == 24)
s_lung_26 <- subset(s_lung, subset = fov == 26)
s_lung_27 <- subset(s_lung, subset = fov == 27)

s_lung_tumor_24 <- subset(s_lung_tumor, subset = fov == 24)
s_lung_tumor_26 <- subset(s_lung_tumor, subset = fov == 26)
s_lung_tumor_27 <- subset(s_lung_tumor, subset = fov == 27)

# 
Idents(s_lung_24) <- "NS_Insitu_Celltype.level1"
ImageDimPlot(s_lung_tumor_24, fov="PSR01", size=5, cols = "glasbey", axes=TRUE, dark.background = FALSE)
ImageDimPlot(s_lung_tumor_27, fov="PSR01", size=5, cols = "glasbey", axes=TRUE, dark.background = FALSE)

# too far apart
#ImageDimPlot(s_lung, fov="PSR01",  cols = "glasbey", axes=TRUE, dark.background = FALSE)

# FOV
# Niche1 (fef*) looks best for Lung 
DefaultBoundary(s_lung@images$PSR01) <- "segmentation"
Idents(s_lung) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"

# Idents(s_lung_24)

s_lung_24 <- subset(s_lung, subset = fov == 24)
s_lung_27 <- subset(s_lung, subset = fov == 27)
p_lung_24_spat_clusters <- ImageDimPlot(s_lung_24, fov="PSR01", size=2, dark.background = FALSE) + # cols = "glasbey", "polychrome"
  ggtitle("Lung FOV 24 - clusters")
ImageDimPlot(s_lung_27, fov="PSR01", size=2, cols = "glasbey",  dark.background = FALSE) + 
  ggtitle("Lung FOV 27 - niche 1")

p_lung_24_spat_niches <- ImageDimPlot(s_lung_24, fov="PSR01", size=2, dark.background = FALSE) + # cols = "glasbey", "polychrome"
  ggtitle("Lung FOV 24 - clusters")
ImageDimPlot(s_lung_27, fov="PSR01", size=2, cols = "glasbey",  dark.background = FALSE) + 
  ggtitle("Lung FOV 27 - niche 1")

Idents(s_lung_24) <- "NS_Insitu_Celltype"

Idents(s_lung_24) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" # niche1
#Idents(s_lung) <- "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"
p_lung_24_spat_niche1 <- ImageDimPlot(s_lung_24, fov="PSR01", size=2, axes=TRUE, dark.background = FALSE) +
  ggtitle("Lung FOV 24 - niche 1 is close")
Idents(s_lung_27) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" # niche1
p_lung_27_spat_niche1 <- ImageDimPlot(s_lung_27, fov="PSR01", size=2, axes=TRUE, dark.background = FALSE) + 
  ggtitle("Lung FOV 27 - niche 1")

# Eval which niche is best for which? view more broadly
Idents(sobj1) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
ImageDimPlot(sobj1, fov="PSR01", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Niche 1")


# this one has an additional niche class mixed in (not more elucidating)
Idents(sobj1) <- "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"
ImageDimPlot(sobj1, fov="PSR01", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Niche 2")
sobj2 <- subset(sobj, subset = Run_Tissue_name == "PSR-02")
Idents(sobj2) <- "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"
ImageDimPlot(sobj2, fov="PSR02", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Niche 2")

# Niche1 notices tissue diff's Niche 2 not helpful here
Idents(sobj2) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
ImageDimPlot(sobj2, fov="PSR02", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Niche 1")
Idents(sobj2) <- "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"
ImageDimPlot(sobj2, fov="PSR02", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Niche 2")

s_lung_tumor_24 <- subset(s_lung_24, subset = path_annot_id > 0)


#################################
# Patient2 plots
#################################
# Workaround for 
# PrepSCTFindMarkers(s_tumor_2)
# Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 87
# Error in PrepSCTFindMarkers(s_tumor_2) : 
#   Multiple UMI assays are used for SCTransform: Nanostring, RNA
# Object contains multiple models with unequal library sizes

# If Subset with both Organoid and Patient samples won't merge or integrate in Seurat
# use this to build a new non-spatial Seurat object for basic DEA
force_merge <= function(s_subset_w_PDO) {
  
  # dump then reload new super-basic Seurat object with combined cell counts
  counts <- as.matrix(GetAssayData(s_subset_w_PDO, layer = "counts"))
  dim(counts)
  
  df_metadata <- s_subset_w_PDO@meta.data
  s_new <- CreateSeuratObject(counts = counts, meta.data = df_metadata, assay = "RNA")
  dim(s_new)
  s_new <- SCTransform(s_new, assay = "RNA", verbose = FALSE)
  
  # Note: these can be highly customized; run these as defaults but can be re-run outside the function too
  s_new <- RunPCA(s_new, assay = "SCT", verbose = FALSE)
  return(s_new)
}


# dump then reload new object with combined cell counts
s_temp <- s_tumor_2
counts <- as.matrix(GetAssayData(s_temp, layer = "counts"))
dim(counts)

df_metadata <- s_temp@meta.data
s_new <- CreateSeuratObject(counts = counts, meta.data = df_metadata, assay = "RNA")
dim(s_new)
s_new <- SCTransform(s_new, assay = "RNA", verbose = FALSE)
s_new <- RunPCA(s_new, assay = "SCT", verbose = FALSE)
ElbowPlot(s_new)
s_new <- RunPCA(s_new, assay = "SCT", verbose = FALSE,  npcs = 8)
Loadings(s_new, reduction = "pca")[1:5, 1:5]
s_new <- FindNeighbors(s_new, assay = "SCT", dims = 1:8, verbose = FALSE)
s_new <- FindClusters(s_new, resolution = 0.5, verbose = FALSE)
s_new <- RunUMAP(s_new, assay = "SCT", dims = 1:8, verbose = FALSE)
# By cluster for cluster naming
mkrs <- FindAllMarkers(s_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
mkrs <- FindAllMarkers(s_new)

Idents(s_new) <- "Sample.ID"
mkrs.Samples <- FindAllMarkers(s_new, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
mkrs.Samples

Idents(s_new) <- "Sample.Label"
DoHeatmap(s_new, features = VariableFeatures(s_new)[1:50])
ggsave(paste0(results_dir, "/heatmap_Pt2_top50_heatmap.png"), last_plot())
FeaturePlot(s_new, features = gentles_adeno, ncol = 2)
RidgePlot(s_new, features = gentles_adeno, ncol = 2)

# volcano plot
# 2_3 vs 2_4

VolcanoPlot(s_new, group.by = "Sample.Label", ident.1 = "2_3", ident.2 = "2_4", min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay="SCT")


sample.de.markers.Pt2tumor <- FindMarkers(s_new, ident.1 = "2_3", ident.2 = "2_4", assay="SCT", recorrect_umi=TRUE )
mkrs

table(s_tumor_2@meta.data$Sample.Label)

Idents(s_new) <- "Sample.ID"
table(s_new@meta.data$Sample.ID)
sample.de.markers.tumor.t1.t2 <- FindMarkers(s_new, ident.1 = "2_3", ident.2 = "2_4")
sample.de.markers.tumor.t1.ORG <- FindMarkers(s_new, ident.1 = "2_3", ident.2 = "ORG_2_3")
dim(sample.de.markers.tumor.t1.t2)
dim(sample.de.markers.tumor.t1.ORG)
# view results
# reorder so gene is first
sample.de.markers.tumor.t1.t2$Gene <- rownames(sample.de.markers.tumor.t1.t2)
sample.de.markers.tumor.t1.ORG$Gene <- rownames(sample.de.markers.tumor.t1.ORG)

colnames <- colnames(sample.de.markers.tumor.t1.t2)
colnames <- c("Gene", colnames[1:5])
sample.de.markers.tumor.t1.t2 <- sample.de.markers.tumor.t1.t2[,colnames]
sample.de.markers.tumor.t1.ORG <- sample.de.markers.tumor.t1.ORG[,colnames]
head(sample.de.markers.tumor.t1.t2)

write_csv(sample.de.markers.tumor.t1.t2, paste0(results_dir, "/Patient2_t1_t2Tumor_DEGs.csv") )
write_csv(sample.de.markers.tumor.t1.ORG, paste0(results_dir, "/Patient2_t1_ORG_Tumor_DEGs.csv") )


# group.by='seurat_cluster.0.2'
DefaultAssay(s_new) <- "RNA"
p_new <- DimPlot(s_new, assay="SCT") + ggtitle("Lung Tumor only - cell clusters")


for (i in 1:4){
  print(i)
  print(slide_names[i])
  print(run_names[i])
  flatfiles_dir <- paste0(rna_root_dir, '/', slide_names[i], '/', run_names[i])
  print(flatfiles_dir)
  print(file.exists(flatfiles_dir))
  list.files(flatfiles_dir)
  
  # Fix polygons file to re-orient global y coords, do only once. 
  #success <- reorient_polygon_global_y(flatfiles_dir, slide_name = slide_names[i])
  
  obj_list[[i]] <- df_exprMat <- read.csv(paste0(flatfiles_dir, slide_name, '_exprMat_file.csv'))
  # LoadNanostring doesn't handle the metadata 
  obj_list[[i]] <- load_meta_load_missed(obj_list[[i]], paste0(flatfiles_dir, "/",slide_names[i],"_metadata_file.csv.gz"))
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
  obj_list[[i]] = remove_sys_control(obj_list[[i]])
  print(paste0("Num Features after removing System Controls (expect 960): ", nrow(obj_list[[i]])))
  
  #obj_list[[i]][["Nanostring"]]@counts <- remove_sys_control(obj_list[[i]], "Nanostring", "counts")
  
  # Negative controls should be subtracted/removed before normalization, but not before running InSituType
  # However, we keep 0 values for Negative controls in the counts matrix after subtracting, and set aside the means for later
  negmeans_list[[i]] <- get_neg_control_means(obj_list[[i]])
  
  #SetAssayData(obj_list[[i]], layer = "counts", new.data = subtract_neg_control(obj_list[[i]], "Nanostring", "counts"))
  obj_list[[i]] <- subtract_neg_control(obj_list[[i]])
  print(paste0("Num Features after accounting for Negative Controls (s/b unchanged): ", nrow(obj_list[[i]])))
  print(paste0("Num Cells before QC: ", ncol(obj_list[[i]])))
  # BASIC QC filter needed to prevent sparsity errors in SCTransform; option for more later 
  obj_list[[i]] <- subset(obj_list[[i]], subset = (nFeature_RNA > 10 & nCount_RNA > 20)) #20 genes/cell. also vs nFeature_Nanostring
  print(paste0("Num Cells after super basic QC: ", ncol(obj_list[[i]])))
  # Remove cells that didn't pass the AtoMx QC
  obj_list[[i]] <- subset(obj_list[[i]], subset = qcCellsFlagged == FALSE)
  print(paste0("Num Cells after filtering on more complete cell QC: ", ncol(obj_list[[i]])))
  
  # Normalizes counts can be done later after QC
  obj_list[[i]] <- SCTransform(obj_list[[i]], assay = "Nanostring", verbose = TRUE)
  obj_list[[i]]
  
}

# Merge Seurat objects using chaining
# Merging allows for group metadata to be added to the Seurat object, and also allows for
# joint dimensional redux and clustering on the underlying RNA expression data if desired. (but we won't)
sobj <- merge(obj_list[[1]], y = c(obj_list[[2]], obj_list[[3]], obj_list[[4]]))

# flatfiles_dir <- paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_name,'/', run_name, '/')
# flatfiles_dir
# # list the files found in flatfiles_dir
# flatfiles <- list.files(flatfiles_dir)
# flatfiles # like slide_name, '_exprMat_file, '_fov_positions_file', 'metatdata_file', '-polygons.csv'

# todo: this is better as data.table, see cosMx example
#df_exprMat <- read.csv(paste0(flatfiles_dir, slide_name, '_exprMat_file.csv'))
#df_fov_pos <- read.csv(paste0(flatfiles_dir, slide_name, '_fov_positions_file.csv'))
#df_metadata <- read.csv(paste0(flatfiles_dir, slide_name, '_metadata_file.csv'))
#df_polygons <- read.csv(paste0(flatfiles_dir, slide_name, '-polygons.csv'))



# ideas: Integration?  merge Nanostring assays, then Integrate? 
s_tumor_2 <- subset(sobj, subset = PatientID == 2 & NS_Insitu_Celltype.level1 == "tumor") 
s_tumor_2 <- SCTransform(s_tumor_2, assay="Nanostring", verbose = TRUE)

Layers(s_tumor_2)
s_tumor_2 <- RunPCA(s_tumor_2, verbose = FALSE, npcs = 20)
s_tumor_2 <- IntegrateLayers(object = s_tumor_2, method = CCAIntegration, 
                             #orig.reduction = "pca",
                             new.reduction = "integrated.cca",
                          verbose = FALSE)

summary(s_tumor_2@meta.data$nCount_RNA)
# re-join layers after integration
s_tumor_2[["Nanostring"]] <- JoinLayers(s_tumor_2[["Nanostring"]])


s_pt2 <- subset(sobj, subset = Run_Tissue_name == "PSR-01" & PatientID == 2)
s_org2 <- subset(sobj, subset = Run_Tissue_name == "PSR-Organoids" & PatientID == 2)
s_tumor_2 <- merge(s_pt2, s_org2)
s_tumor_2 <- SCTransform(s_tumor_2, assay="Nanostring", verbose = TRUE)

#s_tumor_2 <- subset(s_patient2, subset = NS_Insitu_Celltype.level1 == "tumor")
#s_tumor_2 <- subset(s_tumor_2, subset = fov %in% c(1,2,3,4,5,6,7,8,24,25,27))

# Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-

# tumor cells only (by InsituType)
Idents(s_tumor_2) <- "NS_Insitu_Celltype"
ImageDimPlot( s_tumor_2, fov="PSR01", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Tumor FOV 1")
ImageDimPlot( s_tumor_2, fov="PSROrganoids", size=1, cols = "glasbey", axes=TRUE, dark.background = FALSE) +
  ggtitle("Tumor FOV 1")

s_tumor_2 <- PrepSCTFindMarkers(s_tumor_2, assay = "SCT", verbose = TRUE)
#s_tumor_2 <- SCTransform(s_tumor_2, assay="Nanostring", verbose = FALSE)
DefaultAssay(s_tumor_2) <- "SCT"
# 
s_tumor_2 <- FindVariableFeatures(s_tumor_2, selection.method = "vst", nfeatures = 500)
# assume FindVariableFeatures has been run
p1 <- VariableFeaturePlot(s_tumor_2)
top10 <- head(VariableFeatures(s_tumor_2))
p2 <- LabelPoints(plot = p1, points=top10, repel=TRUE)
p2

s_tumor_2 <- RunPCA(s_tumor_2, verbose = FALSE, npcs = 20)

# Assess how PCA did with different npcs, 20, 30, 50 e.g. 
# ElbowPlot to decide best cutoff.  DimPlot colored by possible batch effect variables. e.g. slides should look integrated.
ElbowPlot(s_tumor_2)
s_tumor_2 <- RunPCA(s_tumor_2, verbose = FALSE, npcs = 6)

DimPlot(s_tumor_2, group.by="Sample.Label")

ImageDimPlot(s_tumor_2, fov="PSR01", group.by="NS_Insitu_Celltype.level1", cols = "glasbey", dark.background = FALSE)
ImageDimPlot(s_tumor_2, fov="PSROrganoids", group.by="NS_Insitu_Celltype.level1", size=5, cols = "glasbey", dark.background = FALSE)


ggsave(paste0(results_dir, "ImageDimPlot_Pt2__4Clusters.png"), last_plot())
ImageDimPlot(s_tumor_2, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_pt2_4Clusters.png"), last_plot())

# same plots with markers
#Idents(s_tumor_2) <- 'EPCAM'
#ImageDimPlot(s_tumor_2, fov="PSR01", group.by="EPCAM", size=.5, cols = "glasbey", dark.background = FALSE) # , features = c('EPCAM', 'PTPRC', 'CD31', 'CD10'))

p1 <- ImageFeaturePlot(s_tumor_2, fov="PSR01", features="EPCAM", size=2, dark.background = FALSE)
p2 <- ImageFeaturePlot(s_tumor_2, fov="PSROrganoids", features="EPCAM", size=2, dark.background = FALSE)
ggsave(paste0(results_dir, "ImageFeaturePlot_pt2_EPCAM.png"), p1 + p2)

# something went wrong with plotting too much here
# Warning: No FOV associated with assay 'SCT', using global default FOV
#ImageDimPlot(s_lung_24, molecules="EPCAM", size=5, mols.cols = "red", dark.background = FALSE)

DimPlot(s_tumor_2, reduction = "pca", group.by = "Sample.ID", alpha=0.5, cols='polychrome')
# What features were used in PCA?

DefaultAssay(s_tumor_2) <- "SCT"
s_tumor_2 <- FindNeighbors(s_tumor_2, reduction="pca", dims = 1:6)
s_tumor_2 <- FindClusters(s_tumor_2, resolution=0.3, cluster.name='seurat_cluster.0.3',verbose = TRUE)

s_tumor_2 <- RunUMAP(s_tumor_2, dims = 1:6)
p_pt2_clusters <- DimPlot(s_tumor_2, group.by='seurat_cluster.0.3') + ggtitle("Lung Tumor only - cell clusters")


DimPlot(s_tumor_2, reduction = "umap")
DimPlot(s_tumor_2, group.by="fov", reduction = "umap")

Idents(s_tumor_2) <- "Sample.ID"
s_tumor_2 <- PrepSCTFindMarkers(s_tumor_2, assay="Nanostring")
mkrs <- FindAllMarkers(s_tumor_2)
sample.de.markers.Pt2tumor <- FindMarkers(s_tumor_2, ident.1 = "2_3", ident.2 = "2_4", assay="SCT", recorrect_umi=TRUE )

table(s_tumor_2@meta.data$Sample.Label)

# view results
# reorder so gene is first
sample.de.markers.tumor$Gene <- rownames(sample.de.markers.tumor)
colnames <- colnames(sample.de.markers.tumor)
colnames <- c("Gene", colnames[1:5])
sample.de.markers.tumor <- sample.de.markers.tumor[, c(3, 1, 2)]
head(sample.de.markers.tumor)

write_csv(sample.de.markers.tumor, paste0(results_dir, "LungTumor_DEGs.csv") )

library(gridExtra)

####################################################
#. DESeg2
###################################################
#library(DESeq2)

# get the count
counts <- as.matrix(GetAssayData(s_lung, layer = "counts"))
head(counts)

# Example metadata with a condition column
metadata <- s_lung@meta.data
metadata$Sample.ID <- factor(metadata$Sample.ID)
condition <- factor(metadata$Sample.ID)  # Ensure your condition column is a factor
condition

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

dds <- DESeq(dds)
results <- results(dds)

# Find DE features between Samples
Idents(s_lung) <- "Sample.ID"
table(s_lung@meta.data$Sample.ID)
sample.de.markers <- FindMarkers(s_lung, ident.1 = "1_1", ident.2 = "2_4")
# view results
head(sample.de.markers)
write_csv2(sample.de.markers, paste0(results_dir, "Lung_DEGs.csv"))


# Cell Type specific:
s_lung_tumor <- subset(s_lung, subset = NS_Insitu_Celltype.level1 %in% c("tumor", "epithelial"))
Idents(s_lung_tumor) <- "Sample.ID"
table(s_lung_tumor@meta.data$Sample.ID)
sample.de.markers.tumor <- FindMarkers(s_lung_tumor, ident.1 = "1_1", ident.2 = "2_4")

# view results
# reorder so gene is first
sample.de.markers.tumor$Gene <- rownames(sample.de.markers.tumor)
colnames <- colnames(sample.de.markers.tumor)
colnames <- c("Gene", colnames[1:5])
sample.de.markers.tumor <- sample.de.markers.tumor[, c(3, 1, 2)]
head(sample.de.markers.tumor)

write_csv(sample.de.markers.tumor, paste0(results_dir, "LungTumor_DEGs.csv") )




####################################################
# Explore the markers we expect
###################################################
lung_24_plots <- list() # of plots

DefaultBoundary(s_lung_24@images$PSR01) <- "segmentation"
# markers of interest, by counts level

# mult features generates 3 plots side by side.  Feature plots can't overlap 
ImageFeaturePlot(s_lung_24, fov="PSR01", features = "EPCAM", dark.background = FALSE)
ImageFeaturePlot(s_lung_24, fov="PSR01", features = immune_markers,  dark.background = FALSE)
ImageFeaturePlot(s_lung_24, fov="PSR01", features = gentles_adeno, dark.background = FALSE)

Idents(s_lung_24) <- "NS_Insitu_Celltype.level1"
lung_24_plots[1] <- ImageDimPlot(s_lung_24, fov="PSR01", dark.background = FALSE) + 
  ggtitle("Lung FOV 24 - NS Insitu Type Level1")

Idents(s_lung_24) <- "NS_Insitu_Celltype.level2"
lung_24_plots[2] <- ImageDimPlot(s_lung_24, fov="PSR01", dark.background = FALSE) + 
  ggtitle("Lung FOV 24 - NS Insitu Type Level2")

# Arrange plots in a 9x9 grid
do.call(grid.arrange, c(plots, ncol = 3, nrow = 4))
# Save the 9x9 grid of plots to a PDF file
pdf("3x3_plots.pdf", width = 20, height = 26)
do.call(grid.arrange, c(plots, ncol = 3, nrow = 3))
dev.off()


####################################################
# Explore the markers our dim redux thinks are important
###################################################

# plot umap w cluser numbers
# heatmap of cluster marker expression: what genes ID each cluster. 

# DimPlot, ImageDimPlot for each marker. 







#################################################
# magick for image overlay
################################################
# Add image?

# open the graphics device for displaying magick images/plots, appears as Viewer tab in RStudio
image_graph(width=400, height=400)
dev.off()
# Example usage of image_graph()
#plot_img <- magick::image_graph(width = 800, height = 600, res = 96) # relies on last_plot
# plot_img <- magick::image_graph(width = image_info(overlay_img)$width, 
#                                 height = image_info(overlay_img)$height, res = 96)
# show(plot_img)
# dev.off() # closes up whatever display devide (Viewer, here)


# Try # /Volumes/T7Shield/PSR-GEN-057/rna/PSR01/20240329_210634_S2/CellStatsDir/Morphology2D/20240329_210634_S2_C902_P99_N99_F024.TIF
# cell seg: /Volumes/T7Shield/PSR-GEN-057/rna/PSR01/20240329_210634_S2/CellStatsDir/Segmentation_dd2fb99a-830d-46bc-be98-38f049830423_001/CellOverlay
# CellOverlay_F024.jpg
cell_overlay_path <- "/Volumes/T7Shield/PSR-GEN-057/rna/PSR01/20240329_210634_S2/CellOverlay"
cell_overlay_file_24 <- paste0(cell_overlay_path, "/CellOverlay_F024.jpg")
cell_overlay_file_27 <- paste0(cell_overlay_path, "/CellOverlay_F027.jpg")
file.exists(cell_overlay_file_27)
overlay_seg_27 <- image_read(cell_overlay_file_27)
overlay_seg_24 <- image_read(cell_overlay_file_24)
plot(overlay_seg_24)

# Try also, Morphology2D 5 channel tifs
morph_path  <- "/Volumes/T7Shield/PSR-GEN-057/rna/PSR01/20240329_210634_S2/CellStatsDir/Morphology2D"
morph_file_24 <- paste0(morph_path, "/20240329_210634_S2_C902_P99_N99_F024.TIF")
morph_file_27 <- paste0(morph_path, "/20240329_210634_S2_C902_P99_N99_F027.TIF")
file.exists(morph_file_27)

overlay_img_24 <- image_read(morph_file_24)  %>% # magick package
  image_convert(type = 'TrueColorAlpha') %>%
  # only black is transparent
  image_transparent("black", fuzz = 0.1)

plot(overlay_seg_24)

#overlay_img # this will start flashing through the 5 channels in the Viewer
length(overlay_img_24) # number of channels
overlay_img_24_1 <- overlay_img_24[1] #PanCK 
overlay_img_24_2 <- overlay_img_24[2] # membrane
overlay_img_24_4 <- overlay_img_24[4] # CD45
overlay_img_24_5 <- overlay_img_24[5] # DAPI
p <- ggplot()
#grid.raster(overlay_img_24_1) # this doesn't work. 
#rasterImage(overlay_img_24_1, 0, 0, 1, 1)
plot(as.raster(overlay_img_24_1))

class(p)
p

overlay_img <- image_read(morph_file_27) # magick package
overlay_img # this will start flashing through the 5 channels
length(overlay_img) # number of channels
overlay_img_27_1 <- overlay_img[1] #PanCK 
overlay_img_27_2 <- overlay_img[2] # membrane
overlay_img_27_4 <- overlay_img[4] # CD45
overlay_img_27_5 <- overlay_img[5] # DAPI
overlay_img_27_5

plot_overlay_24_1 <- plot(overlay_img_24_1)
plot_overlay_24_2 <- plot(overlay_img_24_2)
plot_overlay_24_4 <- plot(overlay_img_24_4)
plot_overlay_24_5 <- plot(overlay_img_24_5)
plot_overlay_24_1

# this opens Viewer
#fig <- image_graph(width = 400, height = 400, res = 96)
s_lung_24 <- subset(sobj1, fov == 24)
p_tumor_annot <- ImageDimPlot(s_lung_24, fov="PSR01",size=2, group.by = "path_annot_name", axes=TRUE, dark.background = FALSE) + 
  ggtitle("FOV 24 - annot")
p_tumor_annot



s_lung_27 <- subset(sobj1, fov == 27)
ImageDimPlot(s_lung_27, fov="PSR01",size=2, group.by = "path_annot_name", axes=TRUE, dark.background = TRUE) + 
  ggtitle("Lung FOV 27 - Pathology annot")


p <- plot_channel_as_raster(overlay_img_24_1, color="green")
p


# Cell Overlay file (color, one layer)

fov <- 24
image_plots <- get_plots_for_fov(fov)
# 1: cell seg overlay, 2: panCK, 4: CD45, 5: DAPI
image_plots[2]

# with cowplot
#cowplot::plot_grid(plot_overlay_24_1, plot_overlay_24_2, plot_overlay_24_4, plot_overlay_24_5, nrow = 2, ncol = 2)
p_feature_dim1 <- FeaturePlot(s_lung_24, features=c("Mean.PanCK", "EPCAM", "KRT7"))
p_feature_dim1

# gridExtra
#class(image_plots[[1]])
p_path_annot_psr1
p_lung_dim_loadings
# molecules view zooms out. 

# group.by = "path_annot_name"

p_epithelial_24 <- ImageFeaturePlot(s_lung_24, fov="PSR01", features=epithelial_markers, size=2, axes=FALSE, 
                                dark.background = FALSE) 
p_epithelial_24

p_macros_24 <- ImageFeaturePlot(s_lung_24, fov="PSR01", features=macrophage_lineage, size=2, axes=FALSE, 
                                dark.background = FALSE) 
p_macros_24

p_immune_24 <- ImageFeaturePlot(s_lung_24, fov="PSR01", features=immune_markers, size=2, axes=FALSE, 
                 dark.background = FALSE) 
p_immune_24
p_epithelial_24
  #ggtitle("FOV 24 - tumor related")

# molecules view: bug with including molecules
ImageFeaturePlot(s_lung_24, fov="PSR01", features="EPCAM", molecules="EPCAM", size=2, axes=FALSE, dark.background = FALSE) + 
  ggtitle("FOV 24 - EPCAM")

gridExtra::grid.arrange(image_plots[[1]], image_plots[[2]], image_plots[[4]], 
                        p_lung_clusters,
                        # overlay known tumor cells and spatial cluster labels
                        p_lung_tumor_clusters,
                        p_tumor_annot,
                        p_lung_24_spat_clusters,
                        #p_lung_24_spat_niche1,
                        # add dim plot, ImageDimPlot for Mean.PanCK, EPCAM, KRT7, etc. 
                      
                  
                        nrow = 3, ncol = 3)





# also scale to fit in viewer

img5 <- overlay_img_5
img5 <- image_normalize(img5)
#img2 <- image_modulate(img2, brightness = 200) # if not bright enough 
img5

# How to create composite by coloring channels?
###############
#image_rgb <- image_convert(overlay_img_5, "rgb")
# only blue
# Set red and green channels to zero
#image_blue <- image_fx(image_rgb, expression = "u.r=0; u.g=0; u.b") #except this doesn't work

# must plot something first. then ADD a raster to the current plot
# rasterImage(img2, 0,400,0,400) # error plot.new() has not been called yet
# I had weird results with this.
#g <- rasterGrob(overlay_img, interpolate = TRUE)
# qplot deprecated (use qqplot), but worked. 
#qplot(1, 1, geom = "blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
#ggdraw() + draw_image(overlay_img, x = 0.5, y = 0.5, width = 1, height = 1)



# Example data and plot
p <- ImageFeaturePlot(s_lung_27, fov="PSR01", features=c("EPCAM","PTPRC"),
                      #cols = c("white", "red"),
                      dark.background = FALSE, size=8) + 
  theme_void() + # Remove background and grid
  theme(
    panel.border = element_blank(),  # Remove borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Remove legends
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),   # Remove y-axis text
    plot.margin = unit(c(0,0,0,0), "cm")  # Remove all margins
  ) 
print(p)
ggsave("temp_plot.png", plot = p, device = "tif", width = 4256, height = 4256, units = "px", dpi = 150,
       bg = "transparent")


Idents(s_lung_27) <- "NS_Insitu_Celltype.level1"
Idents(s_lung_27)

# needs a CROP 
ImageDimPlot(s_lung_27, fov="PSR01", alpha=0.3, molecules= c("EPCAM", "PTPRC"))
  
markers.tumor <- FindMarkers(s_lung_27, ident.1 = "tumor")
# titles  
p2 <- ImageDimPlot(s_lung_27, fov = "PSR01", alpha = 0.3, molecules = rownames(markers.tumor)[1:4],
                   nmols = 10000)
p2
p3 <- ImageDimPlot(s_lung_27, fov = "PSR01", alpha = 0.3, molecules = epithelial_markers,
                        nmols = 10000)
p3

# Convert ggplot to an image
class(overlay_img)
image_info(overlay_img)
image_info(overlay_img)$width
# methods: image_scop, image_scale, image_border, image_trim, image_info, image_read, image_write, image_graph, image_composite, image_annotate, image_convert, image_negate, image_blur, image_charcoal, image_despeckle, image_edge, image_emboss, image_enhance, image_equalize, image_implode, image_median, image_modulate, image_oilpaint, image_quantize, image_rotate, image_sample, image_scale, image_segment, image_sharpen, image_shear, image_threshold, image
# image_modulate # to set brightness, saturation, hue, 
# image_resize (img, "300x300") for example
# layers

# how to save ggplot as raster image in magick, use image_graph()
p
# Render the plot directly to a magick-image object at specified resolution
# open graphics viewer
plot_img <- image_graph(width = 4256, height = 4256, res = 150)
print(p)  # Print the plot to the graphics device opened by image_graph
dev.off()  # Close the device to finalize the image
# Save the image
#image_write(plot_img, path = "saved_plot.png", format = "png") # or ggsave
# Save the plot using ggsave to a temporary file

ImageDimPlot(s_lung_27, fov="PSR01", size=2, cols = "glasbey", axes=TRUE, dark.background = FALSE) 


dev.off() # just in case
plot_img <- image_graph(width = 4256, height = 4256, res = 150)
p <- ImageDimPlot(s_lung_27, fov="PSR01", size=10, cols = "glasbey", dark.background = TRUE) +
  #theme_minimal() +  # Use a minimal theme as a base
  theme(
    panel.border = element_blank(),  # Remove borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Remove legends
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),   # Remove y-axis text
    plot.margin = unit(c(0,0,0,0), "cm")  # Remove all margins
  ) 
print(p) # will be huge
ggsave("temp_plot.png", plot = p, device = "tif", width = 4256, height = 4256, units = "px", dpi = 150,
       bg = "transparent")
dev.off()

# Read the image back using magick
plot_img <- image_read("temp_plot.png")

plot_img <- image_read("temp_plot.png") %>%
  image_convert(type = 'TrueColorAlpha') %>%
  image_transparent("white", fuzz = 0) %>%  # assuming the background is white, changes it to transparent
  image_modulate(brightness = 100, saturation = 100, hue = 100) 

image_attributes(plot_img) # shows alpha is activated

  #image_transparent(opacity = 0.5)  # Set alpha to 50%
# what is this trying to do?
#alpha_channel <- image_channel(plot_img, "alpha") |>
#  image_modulate(brightness = 50, saturation = 100, hue = 100) 

image_info(plot_img)
  
print(plot_img)
#plot_img2 <- image_convert(plot_img, "png", colorspace="RGB", depth=8)
#image_info(plot_img2)

# Optionally delete the temporary file if it's no longer needed
#unlink("temp_plot.png")
#border = 100
l_border = 150
r_border = 80
top_border = 120
bottom_border = 120 # this isn't affecting anything / not stretching
hdim = 4256 - (l_border + r_border)
vdim = 4256 - (top_border + bottom_border)
print(paste("hardcoded margin from ggplot as",border, "by visual inspection"))
# make string out of variables
geometry <- paste(hdim, "x", vdim, "+", l_border, "+", top_border)
geometry
img2 <- image_crop(plot_img, geometry) # "3916x3916+170+170")
#img2
#plot(img2) # to see it all in one screen.

img3 <- (image_scale(img2, "1000x1000!")) # has a border we need removed, as well as a legend. 
#plot(img3)
# Make the plot match AtoMx viewer.  
img4 <- image_flip(img3)
img4 <- image_rotate(img4, 270)
#plot(img4)
#plot_img # huge 

# use grid.raster to draw a raster image within the current grid viewport
#grid.raster(overlay_img, interpolate = FALSE, alpha = 0.5)

# meanwhile, get the background image

#print(image_scale(overlay_img_1, "1000x1000"))
#overlay_img2 <- c(image_scale(overlay_img ,"1000x1000"), image_scale(plot_img, "1000x1000"))
#image_flatten(c(image_scale(overlay_img_1 ,"1000x1000"), image_scale(img2, "1000x1000")), 'Add')
# needs a little manual align.

# This works to combine 2 images from the tiff
#composite_img <- image_composite(image_scale(overlay_img_5, "1000x1000"), image_scale(overlay_img_1, "1000x1000"))
plot(overlay_seg)
overlay_seg2 <- image_convert(overlay_seg, "png", colorspace="RGB", depth=8)
image_info(overlay_seg2)


composite_img <- image_composite(image_scale(overlay_seg2, "1000x1000"), # not matte
                                 image_resize(img4, "1000x1000!"), # matte
                                 
          
                                 operator = "Plus") # note: rescale shrunk the dots.

# Composite the plot over the image
#composite_img <- image_composite(overlay_seg, img2, operator = "over") # looks offset and not blended
#composite_img # huge and doesn't line up (dot plot is smaller)

# View the result. this prints a tibble. 
#class(composite_img)
plot(composite_img) # missing a frame. 

#overlay_seg
dev.off()
# save
image_write(composite_img, paste0(results_dir, "fov_27_level1_composite.png", format="png"))


make_magick_composite_img <- function(overlay_seg_img, plot){
  # Assumes 4256 x 4256 images from cosMx
  # Must first save the ggplot to rasterize it. 
  ggsave("temp_plot.png", plot = plot, device = "png", width = 4256, height = 4256, units = "px", dpi = 300,
         bg = "transparent")

  # Read the image back using magick
  plot_img <- image_read("temp_plot.png")
  plot_img2 <- image_convert(plot_img, "png", colorspace="RGB", depth=8) # transparency option?
  #image_info(plot_img2)
  
  # Optionally delete the temporary file if it's no longer needed
  #unlink("temp_plot.png")
  #border = 100
  l_border = 150
  r_border = 80
  top_border = 120
  bottom_border = 120 
  hdim = 4256 - (l_border + r_border)
  vdim = 4256 - (top_border + bottom_border)
  print(paste("hardcoded margin from ggplot as",l_border,"and", top_border, "by visual inspection"))
  # make string out of variables
  geometry <- paste(hdim, "x", vdim, "+", l_border, "+", top_border)
  geometry
  img2 <- image_crop(plot_img, geometry) # "3916x3916+170+170")
  #plot(img2) # to see it all in one screen.
  
  # geometery with "!" losens the aspect ratio being fixed.
  img2 <- (image_scale(img2, "1000x1000!")) # has a border we need removed, as well as a legend. 
  # Make the plot match AtoMx viewer.  
  img2 <- image_flip(img2)
  img2 <- image_rotate(img2, 270)
  #plot(overlay_seg)
  overlay_seg <- image_convert(overlay_seg, "png", colorspace="RGB", depth=8)
  image_info(overlay_seg)
  
  composite_img <- image_composite(image_scale(overlay_seg, "1000x1000"), 
                                   image_resize(img2, "1000x1000!"),
                                   operator = "Plus") # note: rescale shrunk the dots.
  image_info(composite_img) # magick obj can be plotted, annotated, and saved. 
  return(composite_img)
}

overlay_seg <- image_read(cell_overlay_file_24) %>% 
                          image_convert(type = 'TrueColorAlpha') %>%
                            # only black is transparent
                            image_transparent("black", fuzz = 0.1)
my_plot <- ImageDimPlot(s_lung_24, fov="PSR01", size=10, cols = "glasbey", dark.background = TRUE) +
  #theme_minimal() +  # Use a minimal theme as a base
  theme(
    panel.border = element_blank(),  # Remove borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Remove legends
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),   # Remove y-axis text
    plot.margin = unit(c(0,0,0,0), "cm")  # Remove all margins
  ) 

composite2 <- make_magick_composite_img(overlay_seg, my_plot)
plot(composite2)

Idents(s_lung_24) <- "seurat_cluster.0.3"
ImageDimPlot(s_lung_24, fov="PSR01", size=3, cols = "polychrome", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_lung_fov24_9Clusters.png"), last_plot())

Idents(s_lung_24) <- 'NS_Insitu_Celltype.level1'
ImageDimPlot(s_lung_24, fov="PSR01", size=3, alpha=0.5, cols = "polychrome", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_lung_fov24_level2.png"), last_plot())
# FEATURE PLOT, is bg image supported via fov? 
ImageFeaturePlot(s_lung_24, fov="PSR01", size=2,features="EPCAM") 



Idents(s_lung_24) <- 'NS_Insitu_Celltype.level1'
ImageDimPlot(s_lung_24, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_lung_fov24_level1.png"), last_plot())

ImageDimPlot(s_lung_27, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_lung_fov27_4Clusters.png"), last_plot())

# same plots with markers
Idents(s_lung_24) <- 'EPCAM'
ImageDimPlot(s_lung_24, fov="PSR01", size=5, cols = "glasbey", dark.background = FALSE) # , features = c('EPCAM', 'PTPRC', 'CD31', 'CD10'))

ImageFeaturePlot(s_lung_24, features="EPCAM", size=2) # fov warning but worked 
ImageFeaturePlot(s_lung_24, features="EPCAM", size=2, dark.background = FALSE) # fov warning but worked 
ImageFeaturePlot(s_lung_24, fov="PSR01", features="EPCAM", size=2, dark.background = FALSE)
ggsave(paste0(results_dir, "ImageFeaturePlot_lung_fov24_EPCAM.png"), last_plot())




# Try giving it a 'spatial image'.
# FAILS, need image loaded? where?  
DefaultAssay(s_lung_27)
SpatialFeaturePlot(s_lung_24, features = c("CD4"), images="PSR01") # images = name of image

##########################
# CROP and Coords


# sobj@images$PSR01
# 
# subset_immune_pos <- sobj[, which(genes %in% immune_list)]
# # needs FOV
# # many missing coords?
# df_metadata_s <- subset_immune_pos@meta.data
# str(subset_immune_pos)

fov_inherited <- s_lung_24@images$PSR01
fov_inherited@assay # Nanostring

fov_24a <- s_lung_24@images
DefaultAssay(s_lung_24)
# SpatialFeaturePlot(object = s_lung_24, image = "path/to/your/image.png")
SpatialFeaturePlot(s_lung_24, features = c("EPCAM"), images="PSR01")

df_meta_24 <- s_lung_24@meta.data
coords_df <- df_meta_24[, c("Sample.ID", "CenterY_local_px", "CenterX_local_px")]
head(coords_df) 
# sum(is.na(coords_df$CenterY_global_px)) # lots
# sum(!is.na(coords_df$CenterY_global_px))
# 
dim(coords_df)
# Create from scratch; BUT we'll lose the ploygons?
s_lung_24 <- CreateFOV(s_lung_24, path=cell_overlay_file)
# subset_immune_pos@images$PSR01 <- CreateFOV(sobj@images$PSR01, subset_immune_pos)
# 
# Idents(sobj) <- "Sample.ID"
# ImageDimPlot(subset_immune_pos, fov="PSR01",  cols = "polychrome", dark.background = False)

# better, color bg cells grey, positive cells a color









# Try subset - see if InSituType produces same results on subset as whole?  Yes, they do. 
# sobj1 <- subset(sobj, subset = Run_Tissue_name == "PSR-01")
# sup <- run_insitu_type(sobj1, negmean,title_suffix="PSR-01 NS matrix", group.by="Run_Tissue_name", 
#                        insitu_matrix=insitu_matrix, results_dir=results_dir)
# 



# feature_names <- rownames(sobj[["Nanostring"]])
# # which feature_names in list begin with "Neg"?
# 
# feature_names 
# sys_control_names <- feature_names[grep("^SystemControl", feature_names)]
# # Filter feature names that begin with "Neg"
# neg_feature_names <- feature_names[grep("^Neg", feature_names)]
# print(neg_feature_names)

# Get row indices for neg_feature_names
neg_indices <- which(feature_names %in% neg_feature_names)
non_neg_indices <- which(!(feature_names %in% neg_feature_names | feature_names %in% sys_control_names))
# Create sparse matrices for neg_feature_names and others
neg_matrix <- counts_matrix[neg_indices, ]
non_neg_matrix <- counts_matrix[non_neg_indices, ]

# Verify the results
print(dim(neg_matrix))       # Dimensions of the matrix with neg_feature_names
print(dim(non_neg_matrix))   # Dimensions of the matrix without neg_feature_names


# Plot the AtoMx Cell type breakdowns by slide
# "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"  # first one, supervised, about 30 types
# "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments" # unsupervised, 20 clust
# RNA_nbclust_c5e1d3d8.0a17.443a.994b.9a20f297fd2c_1_clusters"   
# df_meta <- sobj@meta.data %>%
#   mutate(Run_Tissue_name = gsub(" ", "_", Run_Tissue_name)) %>%
#   rename(
#     cell_type = RNA_nbclust_c5e1d3d8.0a17.443a.994b.9a20f297fd2c_1_clusters) %>%
#   #filter(Run_Tissue_name %in% c("AzimuthLung6_12", "CPA_LC_Lines", "HCA", "NSCLC", "safeTME")) %>%
#   select(Run_Tissue_name, cell_type)
# 
# df_celltype <- df_meta %>%
#   group_by(Run_Tissue_name, cell_type) %>%
#   summarise(count = n(), .groups = 'drop')
# 
# p <- plot_celltypes(df_celltype, "Run_Tissue_name", "AtoMx Cell Type Breakdowns by Slide")
# ggsave(paste(results_dir, "AtoMx_Cell_Type_Breakdowns_by_Slide.png", sep = ""), plot=p, width = 10, height = 6, dpi = 300)

#tail(sobj@meta.data$RNA_spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_neighbours_plasma.cell) # integers  



cohort_data <- sobj@meta.data[, colnames(sobj@meta.data) %>%
                                grep ("Mean.PanCK|Mean.CD45|PanCK.PT|CD45.PT|Area.um2", ., value = TRUE)]

colnames(cohort_data)
head(cohort_data) # watch out for NAs
sum(is.na(cohort_data$Mean.PanCK))
#is.na(cohort_data) <- cohort_data == 0

# TODO: apply phenotyping for PanCK +/- and CD45 +/- cells
# For now, use continuous values for clustering
hist(cohort_data[cohort_data$Mean.PanCK])
hist(cohort_data[cohort_data$Mean.CD45 > 50 & cohort_data$Mean.CD45 < 2000,"Mean.CD45"])
hist(cohort_data[cohort_data$Mean.PanCK > 50 & cohort_data$Mean.PanCK < 2000,"Mean.PanCK"])
cohort_data$panCK_PT <- ifelse(cohort_data$Mean.PanCK > 500, 1, 0)
cohort_data$CD45_PT <- ifelse(cohort_data$Mean.CD45 > 50, 1, 0)
cohort_data$PT <- paste0(cohort_data$panCK_PT,'/', cohort_data$CD45_PT)
sobj123 <- AddMetaData(sobj1, cohort_data)
head(sobj123@meta.data)














#######################################
# Spatially plots by Genes
################################
# HERE
# Plot the spatially resolved data

marker_list <- c('CD3E', 'CD4', 'CD8A', 'CD20', 'CD68', 'CD163', 'CD11b', 'CD11c', 'CD14', 'CD15', 'CD16', 'CD56', 'CD66b', 'CD45', 'CD19', 'CD20', 'CD21', 'CD23', 'CD27', 'CD30', 'CD31', 'CD34', 'CD38', 'CD45RA', 'CD45RO', 'CD56', 'CD68', 'CD79a', 'CD83', 'CD86', 'CD123', 'CD138', 'CD163', 'CD206', 'CD209', 'CD235a', 'CD244', 'CD279', 'CD303', 'CD304', 'CD326', 'CD335', 'CD337', 'CD338', 'CD339', 'CD340', 'CD341', 'CD342', 'CD343', 'CD344', 'CD345', 'CD346', 'CD347', 'CD348', 'CD349', 'CD350', 'CD351', 'CD352', 'CD353', 'CD354', 'CD355', 'CD356', 'CD357', 'CD358', 'CD359', 'CD360', 'CD361', 'CD362', 'CD363', 'CD364', 'CD365', 'CD366', 'CD367', 'CD368', 'CD369', 'CD370', 'CD371', 'CD372', 'CD373', 'CD374', 'CD375', 'CD376', 'CD377', 'CD378', 'CD379', 'CD380', 'CD381', 'CD382', 'CD383', 'CD384', 'CD385', 'CD386', 'CD387', 'CD388', 'CD389', 'CD390', 'CD391', 'CD392', 'CD393', 'CD394', 'CD395', 'CD396', 'CD397', 'CD398', 'CD399', 'CD400', 'CD401', 'CD402', 'CD403', 'CD404', 'CD405', 'CD406', 'CD407', 'CD408', 'CD409', 'CD410', 'CD411', 'CD412', 'CD413', 'CD414', 'CD415', 'CD416', 'CD417
                 ')
marker_list <- c('CD3E', 'CD4', 'CD8A')
SpatialPlot(sobj, features = marker_list, cols = c("blue", "red", "green"), 
            label = TRUE, label.size = 3, label.color = "black"
            )
# lists
# immune
immune_list <- c('CD3E', 'CD4', 'CD8A', 'CD20', 'CD68', 'CD163', 'CD11b', 'CD11c', 'CD14', 'CD15', 'CD16', 'CD56', 'CD66b', 'CD45', 'CD19', 'CD20', 'CD21', 'CD23', 'CD27', 'CD30', 'CD31', 'CD34', 'CD38', 'CD45RA', 'CD45RO', 'CD56', 'CD68', 'CD79a', 'CD83', 'CD86', 'CD123', 'CD138', 'CD163', 'CD206', 'CD209')
not_immune_list <- c('FN1', 'EPCAM')

genes <- rownames(sobj)
genes


get_marker_pt <- function(marker, sobj, min_count = 1){
  # Given a marker string for a gene found in sobj count matrix, 
  # return a series of 'Pos' or 'Neg' for each cell as (1/0)
  # min_count set the minimum gene count to be considered positive (may vary)
  
  count_mtx <- GetAssayData(sobj, assay = "SCT", layer = "counts")
  #index <- grep(marker, rownames(count_mtx)) # e.g. index of matching row 150 
  
  # make sure markers are found in our metadata
  if (length(setdiff(marker), rownames(count_mtx)) > 0){
    stop("Marker not found in Seurat count matrix")
  }
  
  index <- which(rownames(count_mtx) == marker)
  print(index)
  # how many cells are positive?
  print(length(which(count_mtx[index,] >= min_count))) # e.g. 8421
  # create new column in the metadata for Neg/Pos status of gene
  vec <- rep(0L, length(count_mtx[index, ]))
  vec[which(count_mtx[index, ] >= min_count)] <- 1
  # return df Series with 'Pos' or 'Neg' (1 or 0)  
  return (vec)
}

mkr_vec <- get_marker_pt('CD4', sobj, min_count = 1)
head(mkr_vec)
formatted_pct <- sprintf("%.2f%%", sum(mkr_vec) / length(mkr_vec) * 100)
print(paste(sum(mkr_vec), " positive cells out of ", length(mkr_vec), " is ", formatted_pct ))
                                                 

# Create dataframe to add to metadata of all PTs

df_metadata <- sobj@meta.data
# loop through all markers
for (marker in genes[1:950]){
  print(marker)
  print(head(mkr_vec))
  mkr_vec <- get_marker_pt(marker, sobj, min_count = 1)

  formatted_pct <- sprintf("%.2f%%", sum(mkr_vec) / length(mkr_vec) * 100)
  print(paste(sum(mkr_vec), " positive cells out of ", length(mkr_vec), " is ", formatted_pct ))
  # add to metadata 
  df_metadata[[marker]] <- mkr_vec
}  
head(df_metadata)
sobj <- AddMetaData(sobj, df_metadata)

pt_colors = c("grey", "blue")


## Single marker example:
# Idents(sobj) <- "CD3E"
# p <- ImageDimPlot(sobj, fov="PSR01", split.by= "Sample.Label", size=01.0, 
#                   cols = pt_colors, crop=TRUE, dark.background = FALSE, alpha = 0.5
# ) + ggtitle("CD3E")
# show(p)
# ggsave(paste0(results_dir, "CD3E.png"), plot = last_plot(), width = 10, height = 10, dpi = 300)


# get percent positive for a marker
get_percent_pos <- function(df_meta, marker){
  formatted_pct <- sprintf("%.2f%%", sum(df_meta[[marker]]) / length(df_meta[[marker]]) * 100)
  return(formatted_pct)
}
get_percent_pos(df_metadata, "CD3E")

#########################################################
# Plot each marker spatially, full slide
########################################################

# loop through all markers
#for (slide_name in slide_names)
for (slide_name in c("PSR01"))  {
  for (marker in genes[1:950]){
    print(marker)
    Idents(sobj) <- marker
    p <- ImageDimPlot(sobj, fov=slide_name, split.by="Sample.Label", size=01.0, 
               cols = pt_colors, crop=TRUE, dark.background = FALSE, alpha = 0.5
               ) + ggtitle(paste(marker, " (", get_percent_pos(df_metadata, marker), ")"))
    #show(p)
    ggsave(paste0(results_dir,"/SingleRNAPlots/", slide_name,"_", sanitize_name(marker), ".png"),
           plot = last_plot(), width = 10, height = 10, dpi = 300)
  }
}


########################################################


#. Some groupings:





# How to make low populations more apparent?

slide_name = "PSR01"
feature_name <- "SLPI"

# Plot the background cells
pt_colors = c("grey", "blue")
p_background <- ImageDimPlot(sobj, fov=slide_name, size= 1.0, 
                             cols = pt_colors, dark.background = FALSE, alpha = 0.5) +
  #scale_color_manual(values = c("gray80")) +  # Making background cells gray
  theme_minimal() +
  ggtitle("Tissue Sample with Highlighted Cell Population")
show(p_background)

# Subset or identify the rare cell population
feature_name = 'SLPI'
gene_counts <- GetAssayData(sobj, layer = "counts")[feature_name,,drop=FALSE ] 
#rare_cells_sobj <- subset(sobj, subset = [[feature_name]]  >= 1)
# using [ slice syntax]
col_sums <- colSums(gene_counts)

# Subset the matrix to retain columns where the sum is greater than or equal to 1
filtered_columns <- col_sums >= 1
# this returns a numeric vector instead of a sparse matrix, named with cell_id 
cells <- names(gene_counts[,filtered_columns])

rare_cells_sobj <- sobj[,cells]

# Plot with rare cells highlighted
p_rare_cells <- ImageDimPlot(rare_cells_sobj, fov=slide_name, size=1.0, cols = "red", dark.background = FALSE) +
  theme_void()  # No background/theme to make it overlay-friendly

# Combine the plots

### Why are x and y reversed?  dunno why but flip them.
final_plot <- p_background + 
  geom_point(data = p_rare_cells$data, aes(x = y, y = x, color = "red"), size=0.1, inherit.aes = FALSE)
#show(final_plot)

# Optionally add a legend or further customize
final_plot <- final_plot + 
  labs(color = "Population") +
  #scale_color_manual(values = c("red")) +
  theme(legend.position = "right")

show(final_plot)



#### put in a function
plot_marker_on_greybg <- function(sobj, slide_name, feature_name, pt_colors, rare_threshold = 1) {
  # Given a Seurat object, slide name, marker, and colors, 
  # plot the background cells with the rare cell population highlighted 
  # return the ggplot
  
  p_background <- ImageDimPlot(sobj, fov=slide_name, size= 1.0, 
                               cols = pt_colors, dark.background = FALSE, alpha = 0.5) +
    #scale_color_manual(values = c("gray80")) +  # Making background cells gray
    theme_minimal() +
    ggtitle("Tissue Sample with Highlighted Cell Population")
  show(p_background)
  
  # Subset or identify the rare cell population
  
  gene_counts <- GetAssayData(sobj, layer = "counts")[feature_name,,drop=FALSE ] 
  rare_cells_sobj <- subset(sobj, subset = `feature_name`  >= 1)
  
  # Plot with rare cells highlighted
  p_rare_cells <- ImageDimPlot(rare_cells_sobj, fov=slide_name, size=1.0, cols = "red", dark.background = FALSE) +
    theme_void()  # No background/theme to make it overlay-friendly
  
  # Combine the plots
  
  ### Why are x and y reversed?  dunno why but flip them.
  final_plot <- p_background + 
    geom_point(data = p_rare_cells$data, aes(x = y, y = x, color = "red"), size=0.1, inherit.aes = FALSE)
  #show(final_plot)
  
  # Optionally add a legend or further customize
  final_plot <- final_plot + 
    labs(color = "Population") +
    #scale_color_manual(values = c("red")) +
    theme(legend.position = "right")

  return(final_plot)
  
}

plot_marker_on_greybg(sobj, "PSR01", "SLPI", pt_colors = c("grey", "blue"))
# subset then plot works
s_lung_1_1 <- subset(sobj, subset = Sample.ID == "1_1")
s_lung_2_4 <- subset(sobj, subset = Sample.ID == "2_4")
plot_marker_on_greybg(s_lung_2_4, "PSR01", "SLPI", pt_colors = c("grey", "blue"))

# SAmple FOVs to zoom in. 
# 


# function to remake Sample FOV in a subset of Seurat sobj




#Idents(sobj) <- "Sample.ID"
#SpatialFeaturePlot(sobj, features = c("YES1"))
# Error in data[rownames(x = coordinates), features[j], drop = FALSE] : 
#subscript out of bounds

# Also to SpatialFeaturePlot
# SpatialFeaturePlot(sobj, features = "CD3E", cols = pt_colors, label = TRUE, label.size = 3, label.color = "black")
# SpatialDimPlot
# SpatialDimPlot(sobj, features = "CD3E", cols = pt_colors, label = TRUE, label.size = 3, label.color = "black")


# Given a list of markers, return a vector of Pos/Neg for a cell if *any* are positive (1)     
# get_marker_grp_pt <- function(markers, sobj, label){
#   # Given a list of markers to group and a lable, save the "OR: 
#   df_meta <- sobj@meta.data
# 
#   # logical OR of df_meta[[markers]]
#   df_meta[[label]] <- as.integer(rowSums(df_meta[[markers]] >= 1) > 0)
#   
#   return (df_meta[[label]])
# }
  

# e.g. get "immune" markers

#################################
# Marker lists
################################
# Immune as CD45+/EpCAM-


# Given a list of markers, return a vector of Pos/Neg for a cell if *any* are positive (1)     
get_marker_grp_pt <- function(df_meta, pos_markers, neg_markers){
  # Given a list of metadata with pos/neg, get vector of new phenotype: 
  
  # make sure markers are found in our metadata
  if (length(setdiff(pos_markers, colnames(df_meta))) > 0){
    print(setdiff(pos_markers, colnames(df_meta)))
    stop("One or more positive markers not found in metadata")
  }
  if (length(setdiff(neg_markers, colnames(df_meta))) > 0){
    print(setdiff(neg_markers, colnames(df_meta)))
    stop("One or more negative markers not found in metadata")
  }
  
  vec <- as.integer (
    (rowSums(df_meta[pos_markers] > 0) > 0) & # Any of the pos_markers greater than 0, AND
    (rowSums(df_meta[neg_markers] == 0) == 2) # All of the neg_markers are 0
  )
  return (vec)
}

# starting a new dataframe of columns to add to metadata
df_meta2 <- df_meta['cell_id']
immune_level_1_pos = c('PTPRC')
immune_level_1_neg = c('EPCAM')
df_meta <- sobj@meta.data
immune_vec <- get_marker_grp_pt(df_meta, immune_level_1_pos, immune_level_1_neg)
df_meta2[['Immume_L1']] <- immune_vec

tumor_level_1_pos = c('EPCAM')
tumor_level_1_neg = c('PTPRC')
vec <- get_marker_grp_pt(df_meta, tumor_level_1_pos, tumor_level_1_neg)
df_meta2[['Tumor_L1']] <- vec

# More groupings:
endothelial_l1_pos = c('CD31')
endothelial_l1_neg = c('EPCAM', 'PTPRC')

fibroblast_l1_pos = c('CD10')
fibroblas_l1_neg = c('EpCAM', 'PTPRC', 'CD31')



# count_mtx <- GetAssayData(sobj, assay = "SCT", slot = "counts")
# grep('CD3E', rownames(count_mtx)) # 150 
# # how many cells are positive?
# # only the counts column for the marker
# #count_mtx[150,]
# 
# length(which(count_mtx[150, ] != 0))  # 8421
# # create new column in the metadata for Neg/Pos status of gene
# dat@meta.data$CD3E <- 'Neg'
# dat@meta.data$CD3E[which(count_matrix[150, ] != 0)] <- 'Pos'


# sobj@images$PSR01
# 
# subset_immune_pos <- sobj[, which(genes %in% immune_list)]
# # needs FOV
# # many missing coords?
# df_metadata_s <- subset_immune_pos@meta.data
# str(subset_immune_pos)
# coords_df <- df_metadata_s[, c("Sample.ID", "CenterY_global_px", "CenterX_global_px")]
# head(coords_df) # NA
# sum(is.na(coords_df$CenterY_global_px)) # lots
# sum(!is.na(coords_df$CenterY_global_px))
# 
# dim(coords_df)
# fov_6_5 <- CreateFOV(coords = coords_df, type = "centroids", assay="RNA", key="PSR-02")
# subset_immune_pos@images$PSR01 <- CreateFOV(sobj@images$PSR01, subset_immune_pos)
# 
# Idents(sobj) <- "Sample.ID"
# ImageDimPlot(subset_immune_pos, fov="PSR01",  cols = "polychrome", dark.background = False)

# better, color bg cells grey, positive cells a color



#####################################
# UMAP Project subset
##############################
# Filter the subset
# subset_data <- dataframe %>%
#   filter(group1 == 'fruit') %>%
#   select(-column, -group1)  # Remove non-numeric columns
# # Run UMAP
# umap_result <- umap(subset_data)
# # Extract the UMAP coordinates
# umap_coords <- as.data.frame(umap_result$layout)
# # Combine the UMAP coordinates with the original subset data if needed
# umap_data <- bind_cols(subset_data, umap_coords)
# # Plot the UMAP result
# ggplot(umap_data, aes(x = V1, y = V2)) +
#   geom_point() +
#   labs(title = "UMAP projection of subset data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal()


# check on x and y
colnames(sobj@meta.data)
sum(is.na(sobj@meta.data$CenterX_local_px)) # 0 


# Save the Seurat object 
saveRDS(sobj, file = paste0(results_dir, "/sobj_from_Load_v2.RDS"))

colnames(sobj@meta.data)
table(sobj@meta.data$"spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments") # niche1 - 9

# "spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_assignments"  
Idents(sobj) <- "qcCellsFlagged"

# Why did certain areas have more QC fail? 
Idents(sobj) <- "qcFlagsFOV"
Idents(sobj) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
Idents(sobj) <- "Sample.ID"
ImageDimPlot(sobj, fov="PSR01", cols="polychrome", dark = FALSE, axes=TRUE, flip_xy = FALSE)
ImageDimPlot(sobj, fov="PSR02", cols="polychrome")
ImageDimPlot(sobj, fov="PSRCPA", cols="polychrome")
ImageDimPlot(sobj, fov="PSROrganoids", cols="polychrome")

#cols <- c("polychrome", "nCount_Spatial", "nFeature_Spatial", "percent.mito", "percent.ribo", "percent.hb")
df_cells_failed <- sobj@meta.data[sobj@meta.data$qcCellsFlagged == TRUE, 1:70]
dim(df_cells_failed) # 19K
table(df_cells_failed$qcFlagsCellPropNeg)
table(df_cells_failed$qcFlagsCellCounts) # 18,579 fail (most are this). limit 20
hist(df_cells_failed$nCountPerCell, breaks=1000)

table(df_cells_failed$qcFlagsRNACounts)
table(df_cells_failed$qcFlagsFOV)
table(df_cells_failed$qcFlagsCellComplex) # 5798 fail for this: means each RNA is a distinct gene.
table(df_cells_failed$qcFlagsCellArea)
table(df_cells_failed$qcFlagsCellCounts)

hist(df_cells_failed$nCountPerCell, breaks=1000)
df2 <- df_cells_failed[df_cells_failed$qcFlagsCellCounts == 'Pass', ] 
hist(df2$nCountPerCell, breaks=100)

# what fovs do we have?
sobj@images
# example of each of 4 FOVs. Why was first one's key = 'Nanostring_'
# $PSR01
# Spatial coordinates for 19909 cells and 1157 molecules
# First 10 molecules: AATK, ABL1, ABL2, ACACB, ACE, ACKR1
# Default segmentation boundary: centroids 
# 1 other segmentation boundaries present: segmentation 
# Associated assay: Nanostring 
# Key: psr01_ 

sobj@images$PSR01

#######################################
# Broken ggplot
#######################################

# viridis plot of relative expression
image_dim_plot_counts <- function(obj, fov, features) {
  # Plot the expression values of the features in the specified FOV
  # use viridis color palette to highlight high vs low.
  
  fov_obj <- subset(sobj, subset = (fov == fov & EPCAM > 0))
  dat <- FetchData(fov_obj, vars = c("EPCAM"), layer="counts")
  # need x and y #  "CenterX_local_px", "CenterY_local_px"
  coords_x <- fov_obj@meta.data$CenterX_local_px
  coords_y <- fov_obj@meta.data$CenterY_local_px
  dat$x <- coords_x
  dat$y <- coords_y
  
  # how does ggplot know how to apply color to the EPCAM column?
  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point() +
    scale_color_viridis_d() # Use for discrete variables
    labs(title = "Relative Expression, Spatial for FOV ") +
    theme_minimal()
  
 show(p)
 return(p)
}

p <- image_dim_plot_counts(s_lung_24, fov=24, features=c("EPCAM")) 
show(p)




DefaultAssay(sobj) # Nanostring
Layers(sobj)
Assays(sobj)

Idents(sobj) <- 'Run_Tissue_name'
plot1 <- VlnPlot(sobj, features = "nCountPerCell", pt.size = 0, log = TRUE) + NoLegend()
show(plot1)
# fix... 
sobj$log_nCount_Spatial <- log(sobj$nCountPerCell)
plot2 <- SpatialFeaturePlot(sobj, features = "log_nCount_Spatial") + theme(legend.position = "right") # error
wrap_plots(plot1, plot2)

# BUG: Seurat Spatial Error #2695
# also check if rownames(metadata) == Cells(object)
all(rownames(sobj@meta.data) == Cells(sobj)) # TRUE

SpatialPlot(sobj, features="AATK") # error
Features(sobj)
# try ImageFeaturePlot instead
ImageFeaturePlot(sobj, features=c("VIM"), fov="PSR01")
# KRTs. This makes a really ugly plot btw.
krts <- c("KRT1", "KRT10" ,"KRT13","KRT14", "KRT15",    
          "KRT16" ,"KRT17","KRT18", "KRT19" ,"KRT20",      
          "KRT23","KRT4","KRT5" ,"KRT7"   ,        # "KRT6A/B/C"
          "KRT8", "KRT80" ,"KRT86" )
ImageFeaturePlot(sobj, features=krts, fov="PSR01")

ImageFeaturePlot(sobj_spleen1, features=c("Cd3d","Cd3e","Cd3g"), fov="fov_31")
ImageFeaturePlot(sobj, features=c("CD4","CD8A","CD8B"), fov="PSR01")


# flatfiles_dir <- paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_name,'/', run_name, '/')
# flatfiles_dir
# # list the files found in flatfiles_dir
# flatfiles <- list.files(flatfiles_dir)
# flatfiles # like slide_name, '_exprMat_file, '_fov_positions_file', 'metatdata_file', '-polygons.csv'

# todo: this is better as data.table, see cosMx example
#df_exprMat <- read.csv(paste0(flatfiles_dir, slide_name, '_exprMat_file.csv'))
#df_fov_pos <- read.csv(paste0(flatfiles_dir, slide_name, '_fov_positions_file.csv'))
#df_metadata <- read.csv(paste0(flatfiles_dir, slide_name, '_metadata_file.csv'))
#df_polygons <- read.csv(paste0(flatfiles_dir, slide_name, '-polygons.csv'))

# flatfiles_dir1 <- paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_name,'/', run_name, '/')
# df_metadata1 <- read.csv(paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_names[1],'/', slide_names[1], '_metadata_file.csv.gz'))
# df_metadata2 <- read.csv(paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_names[2],'/', slide_names[2], '_metadata_file.csv.gz'))
# df_metadata3 <- read.csv(paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_names[3],'/', slide_names[3], '_metadata_file.csv.gz'))
# df_metadata4 <- read.csv(paste0('/Volumes/T7Shield/',exp_name, '/rna/', slide_names[4],'/', run_name, '/', slide_names[4], '_metadata_file.csv'))

# I should have the metadata flat files loaded
#plot qc fields...
plot(qcFlagsRNACounts)
plot(qcCllesFlagged)
hist(nCountPerCell)








sobj <- AddMetaData(sobj, df_metadata2) 
saveRDS(sobj, file = paste0(results_dir, "/sobj_from_Load_v3.RDS"))

#############################################
# Plots using FOV and sample meta
########  
# START HERE!

readRDS()
# add samples as FOVs
# cell typing cont.. 


Idents(sobj) <- 'Sample.ID'
p <- DimPlot(sobj, reduction="umap", group.by= c("tissue", "Run_Tissue_name"))
show(p)
ImageDimPlot(sobj, fov="PSR02" ) # group.by="Run_Tissue_name") # looks terrible




###############################################
# Initial Split up our slides for separate cell typing
##############################################

sobj_CPA <- subset(sobj, subset = Run_Tissue_name == 'PSR-CPA')

colnames(sobj_CPA@meta.data)
# plot clusters spatially
Idents(sobj_CPA) <- CellType_full15
DimPlot(sobj_CPA, reduction="umap", group.by = "CellType_full15")


# set colors
colors_tumor1
cpa_colors = colors_tumor1
names(cpa_colors)
names(cpa_colors) <- c("i", "h", "d", "o", "g", "b", "Plastmacytoid.dendritic.cell", "j", 
                       "Tumor1a", "Tumor1",  "Tumor2" ,
                       "Tumor3",  "Tumor4",  "Tumor5" , "Tumor6" , "Tumor7",  "Tumor8" )
# append missing ones

ImageDimPlot(sobj_CPA, fov="PSRCPA", group.by="CellType_full15", cols=cpa_colors, dark.background = FALSE)

# in sobj_CPA, rename values in CellType_full15 from i to "A1", h to "A2", etc
df_meta_CPA <- sobj_CPA@meta.data

# IF count <10, remove (call it 'other')
# which labels have count < 10
# Get the values in CellType_Labeled that have a count < 10
low_count_values <- df_meta_CPA %>%
  group_by(CellType_Labeled) %>%
  summarise(count = n()) %>%
  filter(count < 10) %>%
  pull(CellType_Labeled)

# Print the list of values
print(low_count_values)


df_meta_CPA <- df_meta_CPA %>%
  mutate(CellType_Labeled = case_when(
    CellType_full15 == 'i' ~ 'A1',  # Set 'i' to 'A1'
    CellType_full15 == 'h' ~ 'A2',  
    CellType_full15 == 'd' ~ 'A3',
    CellType_full15 == 'o' ~ 'B1',
    CellType_full15 == 'g' ~ 'B2',
    CellType_full15 == 'mix' ~ 'B3',
    #CellType_full15 == 'g and f?' ~ 'C1',
    CellType_full15 == 'f' ~ 'C1',
    CellType_full15 == 'j' ~ 'C2',
    CellType_full15 == "Mast.cell" ~ "Other",
    CellType_full15 == "NK.cell" ~ "Other",  
    CellType_full15 == "Neutrophil" ~ "Other",         
    CellType_full15 == "Plasmacytoid.dendritic.cell" ~ "Other",
    CellType_full15 == "T.cell.CD8"  ~ "Other",               
    CellType_full15 == "n"  ~ "Other",
    TRUE ~ CellType_full15  # Keep existing values for other cases
  ))

table(df_meta_CPA$CellType_Labeled)
df_meta2 <- as.data.frame(df_meta_CPA[,c("cell_id","CellType_Labeled")])
sobj_CPA <- AddMetaData(sobj_CPA, metadata = df_meta2)

Idents(sobj_CPA) <- 'CellType_Labeled'
ImageDimPlot(sobj_CPA, fov="PSRCPA", group.by="CellType_Labeled", cols="glasbey")
# label cluster 15s:
# A1 = i, A2 = h, A3 = d, B1 = o, B2 = g, B3 = mix, C1 = g and f?, C2 = j 


##########################
# QC Cleanup, only Patient Cells
##########################

# show rows where remove_flagged_cells is TRUE
#nrow(df_metadata[df_metadata$remove_flagged_cells == TRUE,]) # eg 80
hist(df_metadata$nCount_RNA, breaks = 100)
hist(df_metadata$nFeature_RNA, breaks = 100)
table(df_metadata$qcCellsFlagged) # eg 80

# what are some failing FOVs? 
table(df_metadata$qcFlagsFOV)


#max(df_metadata_f002$cell_ID) # c_1_11_170
# max cell id is c_3_2_99.  how does this work

table(df_metadata$Run_Tissue_name)

head(df_metadata[df_metadata$fov == 10, 2]) # c_2_1_1, c_2_1_2, c_2_1_3, c_2_1_4, c_2_1_5, c_2_1_6)

# bar plot of cells grouped by Run_Tissue_Name and fov 
ggplot(df_metadata, aes(x = Run_Tissue_name, fill = df_metadata$fov)) + geom_bar()
ggplot(df_metadata, aes(x = fov, fill = qcCellsFlagged)) + geom_bar()
ggplot(df_metadata, aes(x = fov, fill = qcFlagsFOV)) + geom_bar()
table(df_metadata$fov)

# do qcFlagsFOV fail a fov completely?
table(df_metadata$fov, df_metadata$qcFlagsFOV)
table(df_metadata$Run_Tissue_name, df_metadata$fov)

# distinct values for df_metadata$fov
unique(df_metadata$fov) # 1, 2, 3, 4, 5, 6, 7, 8, 9, .. 350
# show counts by df_metadata$fov
table(df_metadata$fov) # 1: 1541, 2: 1541, 3: 1541, 4: 1541, 5: 1541, ..

# Method 2: Using summary
print(summary(df_metadata$qcCellsFlagged))

# 
# # Method 3: Using dplyr
# factor_summary <- df_metadata %>%
#   group_by(fov, qcCellsFlagged) %>%
#   summarize(count = n(), .groups="keep") %>%
#   arrange(fov, qcCellsFlagged)
# print(factor_summary)
# 
# # this filters the df but not seurat
# df_qc <- df_metadata[df_metadata$qcCellsFlagged == FALSE,]
# filtered_cells <- df_qc$cell_ID
# 
# # subset the sobject
# sobj <- subset(sobj, cells = filtered_cells)
# sobj@version # 4.1.3

# Save the Seurat object 
saveRDS(sobj, file = paste0(results_dir, "/seurat_v5_object.rds"))
# Optionally, reload the saved object to confirm the upgrade (no)
#sobj <- readRDS(paste0(results_dir,"/seurat_v5_object.rds"))
sobj@version 


Idents(sobj_CPA) <- 'Sample.ID'
colnames(sobj_CPA@meta.data)
head(sobj_CPA@meta.data$CenterX_global_px)
# Rename spatial coordinate columns to standard names
#sobj_CPA <- RenameCells(sobj_CPA, new.names = c("imagerow" = "CenterY_global_px", "imagecol" = "CenterX_global_px"))


###########################
#. Create FOV
# when we subset a seurat object, check first if the FOV is still there.  subsetting to a Sample will keep the FOV. 
###########################

create_fov <- function(sobj, x_col, y_col, molecules, assay="RNA", name="name", subset = NULL) {
  # Assumes the Seurat object has metadata columns like CenterX_global_px and CenterY_global_px
  # creates an FOV for all the cell locations in the Seurat object
  sobj@meta.data$imagerow <- sobj@meta.data[x_col]
  sobj@meta.data$imagecol <- sobj@meta.data[y_col]
  # Ensure that spatial coordinates are set correctly
  
  df_coords <- sobj@meta.data[, c("imagerow", "imagecol")]
  
  # if (!is.null(subset)) {
  #   sobj <- subset(sobj, cells = subset)
  # }
  
  fov <- CreateFOV(coords = df_coords, molecules = molecules, type = "centroids", assay="SCT", name=name)
  sobj@images$PSR01_B <- fov
  #sobj@images <- fov
  #DefaultAssay(sobj@images$image) <- assay
  return(sobj)
}

sobj_CPA <- create_fov(sobj_CPA,  x_col = "CenterX_global_px", y_col = "CenterY_global_px")
# check on it
sobj_CPA@images$image

sobj123 <- create_fov(sobj123, x_col = "CenterX_global_px", y_col = "CenterY_global_px")

# this doesn't work.
SpatialFeaturePlot(sobj_spleen1, crop=TRUE,
                    features = c("Cd3d")) # , cols = c("blue", "red"), label = TRUE)

SpatialDimPlot(sobj123)


sobj_CPA@meta.data$imagerow <- sobj_CPA@meta.data$CenterY_global_px
sobj_CPA@meta.data$imagecol <- sobj_CPA@meta.data$CenterX_global_px
# Ensure that spatial coordinates are set correctly
sobj_CPA@images$image$coordinates <- sobj_CPA@meta.data[, c("imagerow", "imagecol")]


fov_cpa <- CreateFOV(coords = sobj_CPA@images$image$coordinates, type = "centroids", assay="RNA", name="CPA")
sobj_CPA@images$image <- fov_cpa
DefaultAssay(sobj_CPA@images$image) <- "RNA"

Layers(sobj_CPA)
Assays(sobj_CPA)
DefaultAssay(sobj_CPA) <- "RNA"
#convert_v3_to_v5(sobj_CPA)
#DefaultAssay(sobj_CPA) <- "RNA_normalized_cacc3861.a7b3.47a9.90f1.795a3d4cdcad_1"
p <- ImageDimPlot(sobj_CPA, axes = TRUE, cols = "glasbey") 
# No compatible spatial coordinates present
show(p)

head(sobj_CPA@meta.data$Line)






#################################
# Generate CPA Cell Matrix
################################

feature_names <- rownames(sobj_CPA[["RNA"]])
counts_matrix <- sobj_CPA[["RNA"]]@counts
# which feature_names in list begin with "Neg"?

feature_names 
# Filter feature names that begin with "Neg"
neg_feature_names <- feature_names[grep("^Neg", feature_names)]
print(neg_feature_names)

# Get row indices for neg_feature_names
neg_indices <- which(feature_names %in% neg_feature_names)
non_neg_indices <- which(!(feature_names %in% neg_feature_names))
# Create sparse matrices for neg_feature_names and others
neg_matrix <- counts_matrix[neg_indices, ]
non_neg_matrix <- counts_matrix[non_neg_indices, ]

# Verify the results
print(dim(neg_matrix))       # Dimensions of the matrix with neg_feature_names
print(dim(non_neg_matrix))   # Dimensions of the matrix without neg_feature_names


# Function to calculate the average by group
calculate_group_averages <- function(matrix, groups) {
  unique_groups <- unique(groups)
  group_averages <- sapply(unique_groups, function(group) {
    group_indices <- which(groups == group)
    group_matrix <- matrix[, group_indices, drop = FALSE]
    rowMeans(group_matrix)
  })
  return(group_averages)
}

cell_type_list <- sobj_CPA$Line
group_averages <- calculate_group_averages(non_neg_matrix, cell_type_list)
group_averages_matrix <- as.matrix(group_averages)

write.csv(as.data.frame(group_averages_matrix), 
          "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_CPA_LC_Lines.csv",
          row.names = TRUE)

# group_averages_matrix <- as.matrix(
#   read.csv("/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_CPA_LC_Lines.csv",
#   row.names = 1))

# create new summary dataframe with counts by fov
# df_fov_counts <- data.frame(table(df_metadata$fov))
# colnames(df_fov_counts) <- c('fov', 'cell_count')
# head(df_fov_counts)
# # add column with the count of rows with remove_flagged_cells == TRUE by fov
# df_fov_counts$flagged_cell_count <- NA
# for (i in 1:nrow(df_fov_counts)) {
#   fov <- df_fov_counts$fov[i]
#   df_fov_counts$flagged_cell_count[i] <- nrow(df_metadata[df_metadata$fov == fov & df_metadata$remove_flagged_cells == TRUE,])
# }
# head(df_fov_counts)


# Add a columnn for good_cells = (cell_count - flagged_cell_count )
# df_fov_counts$good_cell_count <- df_fov_counts$cell_count - df_fov_counts$flagged_cell_count
# df_fov_counts

# todo: maybe drop FOVs with less than 10 good cells.
# we might re-evaluate what makes a good cell..

################################################################
# Phenotypes expected by CoPilot:
# get phenotypes from seurat object
# phenotypes are in the form of a matrix with columns for each cell type and rows for each cell
# the values are the probability that the cell is of that type
# the rownames are the cell_IDs
# the column names are the cell types


# vs what I have: celltypes as a classification.


# In case of failure on AddMetaData, check if the cell names in metadata match the cell names in the Seurat object
# # Check if the cell names in metadata match the cell names in the Seurat object
# if (!all(rownames(df_metadata2) == colnames(sobj))) {
#   # Reorder metadata to match the cell names in the Seurat object
#   df_metadata2 <- df_metadata2[match(colnames(sobj), rownames(df_metadata2)),]
# }
# check rownames
# all(rownames(df_metadata2) == colnames(sobj)) # FALSE (problem)
# colnames(sobj) # e.g. "c_1_122_675" 236796 entries
# rownames(df_metadata2) # e.g. "c_1_122_675" 236322 entries (aha)


#################################################
# Back to the Seurat workflow
#################################################
#  Filter to remove flagged cells



# TODO: filter and render panCK positive cells only, by metadata Mean.panCK > 0.5
# Render various markers whole slide & per sample. 
# normal Seurat workflow

# simplify colnames that are super long, e.g. remove '61b859e8.677d.4cad.9cf5.44e03fe9960e'
# we have many colnames in metadata like ""RNA_spatialclust_61b859e8.677d.4cad.9cf5.44e03fe9960e_1_neighbours_activated.dendritic.cell"  



#################################################
# Run some Seurat steps on Patient Data
#################################################
sobj123@version
DefaultAssay(sobj123)

# Extract the raw counts and normalized data
raw_counts <- GetAssayData(sobj123, slot = "counts", layer = "RNA")
normalized_data <- GetAssayData(sobj123, slot = "data", layer = "RNA")
scale_data <- GetAssayData(sobj123, slot="scale.data", layer="RNA")

# Check if the normalized data slot is filled
if (all(normalized_data == 0)) {
  print("The Seurat object has not been normalized.")
} else {
  print("The Seurat object has been normalized.")
}

# Compare raw counts and normalized data
if (!all(raw_counts == normalized_data)) {
  print("The counts and normalized data differ, indicating normalization has been performed.")
} else {
  print("The counts and normalized data are identical, indicating normalization has not been performed.")
}

# Check if the scale data slot is filled
if (all(scale_data == 0)) {
  print("The Seurat object has not been scaled")
} else {
  print("The Seurat object has been scaled")
}

# Normalize the data
#sobj <- NormalizeData(sobj)
# Using SCTransform which is more appropriate for more complex datasets and where Integration might
#  be desired down the line.
# except this has been done already.
#sobj123 <- SCTransform(sobj123, assay = "RNA")

# These are not done if SCTransform is used
# sobj <- FindVariableFeatures(sobj)
# Scale the data
#sobj <- ScaleData(sobj)

sobj123 <- RunPCA(object = sobj123, npcs=15)
sobj123 <- FindNeighbors(object = sobj123, dims = 1:10)
sobj123 <- FindClusters(object = sobj123, resolution=0.3, cluster.name='seurat_cluster',)
sobj123 <- RunUMAP(object = sobj123, dims = 1:5)
DimPlot(object = sobj123, reduction = "umap")

# Did we select an appropriate number of dimenstions for PCA?
# Extract standard deviations of principal components
pca_results <- sobj123[["pca"]]
stdev <- pca_results@stdev

# Calculate variance explained
variance_explained <- stdev^2 / sum(stdev^2)

# Create a data frame for plotting
scree_data <- data.frame(
  PC = 1:length(variance_explained),
  VarianceExplained = variance_explained
)
sum(variance_explained[1:50])
sum(variance_explained[1:10])

# Generate the scree plot
scree_plot <- ggplot(scree_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance Explained") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the scree plot
print(scree_plot)


# Check if PCA has been run
if ("pca" %in% names(sobj@reductions)) {
  print("Regular PCA has already been run.")
} else {
  print("Regular PCA has not been run yet.")
}

# Scree Plot? Elbow?

#sobj <- RunPCA(sobj, 30)

# Check if any reduction names start with 'approximatepca'
reduction_names <- names(sobj123@reductions)
approx_pca_exists <- any(grepl("^approximatepca", reduction_names))
approx_pca_exists
reduction_names[grepl("^approximatepca", reduction_names)]

if ("approximatepca" %in% names(sobj123@reductions)) {
  print("Approximate PCA has already been run in AtoMx.")
} else {
  print("Approximate PCA has not been run yet.")
}

# but an approximate PCA has been 
approx_pca <- sobj@reductions$approximatepca_e5a7d9b1.8e92.4775.b01c.edc0abb1ae2a_1

# View the reduction object summary
print(approx_pca)

# Inspect embeddings
embeddings <- sobj@reductions$pca@cell.embeddings
print("Embeddings dimensions:")
print(dim(embeddings))
print("First few embeddings:")
print(head(embeddings))

# Inspect feature loadings
loadings <- sobj@reductions$pca@feature.loadings
print("Loadings dimensions:")
print(dim(loadings))
print("First few loadings:")
print(head(loadings))

# Extract PCA loadings
pca_loadings <- sobj123[["pca"]]@feature.loadings

# Extract loadings for PC1
pc1_loadings <- pca_loadings[, 1]

# Create a data frame for plotting
pc1_loadings_df <- data.frame(
  Gene = names(pc1_loadings),
  Loading = pc1_loadings
)

# Sort by absolute value of loadings and select top 20
top20_loadings <- pc1_loadings_df[order(abs(pc1_loadings_df$Loading), decreasing = TRUE), ][1:20, ]

# Plot PCA loadings for PC1 top and bottom 20 
ggplot(top20_loadings, aes(x = reorder(Gene, -Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "PCA Feature Loadings for PC1",
       x = "Gene",
       y = "Loading") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Loadings(sobj123, reduction="pca")

# Perform UMAP
#sobj <- RunUMAP(sobj, dims = 1:10)
# Plot UMAP
DimPlot(sobj123, reduction = "umap")


# Perform UMAP
# Normalize the data
# sobj <- NormalizeData(sobj)
# # Find variable features
# sobj <- FindVariableFeatures(sobj)
# # Scale the data
# sobj <- ScaleData(sobj)
# # Perform PCA
# sobj <- RunPCA(sobj)
# sobj <- RunUMAP(sobj, dims = 1:30)
# # Plot UMAP
# DimPlot(sobj, reduction = "umap")



###########################
# CPA Seurat basics
##########################
sobj_CPA <- SCTransform(object = sobj_CPA, assay="Nanostring")
sobj_CPA <- RunPCA(object = sobj_CPA, npcs=30)
sobj_CPA <- FindNeighbors(object = sobj_CPA, dims = 1:30)
sobj_CPA <- FindClusters(object = sobj_CPA, resolution=0.3, cluster.name='seurat_cluster')
sobj_CPA <- RunUMAP(object = sobj_CPA, dims = 1:10)
p <- DimPlot(object = sobj_CPA, reduction = "umap")
show(p)
ggsave("Umap_CPA_24.png", plot = p)



# plot whole slide with PanCK cells only
# plot whole slide with CD45 cells only
# TODO: how to do spatial plot from Seurat?

df_metadata <- sobj@meta.data
#df_metadata_f002 <- df_metadata[df_metadata$fov == 2,]
#df_metadata_f016 <- df_metadata[df_metadata$fov == 16,]
#filtered_cells <- df_metadata_f016$cell_ID
#sobj_fov16 <- subset(sobj, cells = filtered_cells)

# sobj_fov16 %>%
#   ggplot() +
#   geom_point(aes(x=x_FOV_px, y=y_FOV_px, color = panCK_PT)) +
#   coord_fixed() +
#   labs(titel = "PanCK cells only") +
#   alex_theme()

# have color for each combo PT (panCK, CD45)

col_list <- c("darkblue", "green", "black", "grey", "magenta")

sobj_fov16@meta.data$PT <- paste0(ifelse(sobj_fov16@meta.data$panCK_PT == 1, "PanCK+", "PanCK-"), "/",
                                  ifelse(sobj_fov16@meta.data$CD45_PT == 1, "CD45+", "CD45-"))
cols <- col_list[seq_along(unique(sobj_fov16$PT))]
names(cols) <- unique(sobj_fov16@meta.data$PT)
cols

plot(sobj_fov16$x_FOV_px, sobj_fov16$y_FOV_px, pch = 16, cex = .75, asp = 1, cex.main = 0.75,
     main = "cells in physical space",
     col = cols[sobj_fov16$PT],

)
legend("bottomleft", pch=16, col=cols, legend=names(cols), cex=0.7)
cols
# note: darkblue is both PTs (reduce)


###################################
# Plot by each sample
##################################
subset_1_1 <- subset(sobj, subset = Sample.ID == "1_1")
subset@images$PSR01
ImageDimPlot(subset_1_1, alpha=0.5, cols='polychrome',  dark.background=FALSE)

DimPlot(subset_1_1,  alpha=0.5, cols='polychrome')
#ggsave(paste0(results_dir, "PCA_plot_Run_Tissue_Name123_20pcas.png"), plot = last_plot(), width = 10, height = 10, dpi = 300)


#############################################
# Try plotting whole liver sample 3_6
########################################
df_sample_3_6 <- df_metadata[df_metadata$Sample.ID == "3_6",]
filtered_cells <- df_sample_3_6$cell_ID
sobj_3_6 <- subset(sobj, cells = filtered_cells)


# HERE
# TODO: mult FOVs.
# how to plot subset gracefully, get bounding box of global px coords and subtract offset
plot(sobj_3_6$x_FOV_px, sobj_3_6$y_FOV_px, pch = 16, cex = .75, asp = 1, cex.main = 0.75,
     main = "sample cells in physical space",
    # col = cols[sobj_3_6$PT],
)


                              
cohort_data[1:5, ]
# perform automatic cohorting
cohort <- fastCohorting(cohort_data[,c("Mean.PanCK", "Mean.CD45")],
                        gaussian_transform = TRUE)
# cohort <- fastCohorting(cohort_data[,c("panCK_PT", "CD45_PT")], 
#                         gaussian_transform = TRUE)
table(cohort_data$panCK_PT)
# ("Gaussian_transform = TRUE" maps variables to gaussians in order to 
#  place dramatically different variables on the same scale.)
table(cohort) # 3 cohorts shows possibly good separation

# counts <- sobj123@assays$SCT@counts %>%
#   as.matrix() %>%
#   t()
# 
# counts <- GetAssayData(sobj123, slot = "data", layer="SCT") %>%
#   as.matrix() %>%
#   t() 

dim(counts)
length(negmean)
# error: caused by a rowsum being NA.  avoid this. 
any(rowSums(counts) == 0) # NA values in counts are a problem 
sup <- insitutypeML(x = counts,
                    neg = negmean,
                    cohort = NULL, # cohort_data,
                    reference_profiles = insitu_matrix)   
# 86 genes in the count data are missing from reference_profiles and will be omitted from cell typing (ok-ish)
# Error: Not compatible with requested type: [type=list; target=double]

# ImmuneTumor_safeTME_profilematrix.csv has 693 genes not matching  
# try Lung_HCA_profilematrix

# TODO: provide alias table and rename reference 
# find 'VEGF' in ioprofiles
#ioprofiles['VEGF',]

sup$clus[1:10]
length(sup$clus)

# nice plots
heatmap(sweep(sup$profiles, 1, pmax(apply(sup$profiles, 1, max), .2), "/"), scale="none",
                   main = "/InSituType Clustering - NSCLC public data",
                  )
# heatmap isn't a ggplot so this doesn't work; Save the plot as a file instead.
ggsave(paste0(results_dir,"/ImageDimPlot_NSCLC_heatmap_allslides.png"), plot = p)

iocolors # list of hex codes
cols <- iocolors[seq_along(unique(sup$clust))]
names(cols) <- unique(sup$clust)




fp <- flightpath_plot(flightpath_result = NULL,
                      insitutype_result = sup,
                      col = cols[sup$clust])
print(fp)
ggsave(paste0(results_dir,"ImageDimPlot_NSCLC_FlightPath_all.png"), plot = fp)

length(sup$clust)
sobj123 <- AddMetaData(sobj123, metadata = data.frame(CellType_ns = sup$clust))

Idents(sobj123) <- "CellType_ns"
#Idents(sobj) <- "spatialclust_fef9b11c.70b7.4905.b96d.0b80c5560fb7_1_assignments"
p <- ImageDimPlot(sobj123, fov="PSR01", cols="polychrome", dark.background = FALSE)
show(p)
ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_PSR01.png"), plot = p)
p <- ImageDimPlot(sobj123, fov="PSR02", cols="polychrome", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_PSR02.png"), plot = p)
p <-ImageDimPlot(sobj123, fov="PSROrganoids", cols="polychrome", dark.background = FALSE)
ggsave(paste0(results_dir, "ImageDimPlot_celltype_ns_Organoids.png"), plot = p)



# cell counts by celltype as bar plot
barplot(table(sup$clust), col=cols, las=2, cex.names=0.75, main="Cell counts by InSituType")
# by fov

df_meta <- sobj123@meta.data
df_grp  <- df_meta %>%
  group_by(CellType_safeCPA, Run_Tissue_name) %>%
  summarise(count = n())

df_meta$count <- 1
# Create the barplot
ggplot(df_grp, aes(x = Run_Tissue_name, y=count, fill = CellType_safeCPA)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Counts by Celltype Grouped by Slide",
       x = "Celltype",
       y = "Count") +
  theme_minimal()


# Projected 
RunPCA(sobj1,)
RunUMAP(sobj1, dims = 1:10, seed.use = 123)

plot(sobj123@reductions$umap@cell.embeddings, pch=16, cex=0.75,
     asp=1, cex.main=0.75, 
     main = "cells in UMAP space",
     col=cols[sup$clust],
     xlab="", ylab="", xast="n", yast="n")
legend("bottomleft", pch=16, col=cols, legend=names(cols), cex=0.7)

show(up)
ggsave("ImageDimPlot_Azimuth_umap_PSR01.png", plot = up)


######################################
# Single Sample plot
#####################################


# create FOV











############################
# InSituType Cell Typing for Cell Pellet Array
############################

#insitu_matrix <- as.matrix(read.csv( "/Users/annmstrange/Documents/cosMx/InSituType/CellProfile_matrix_AzimuthLung6_12.csv",
#                                     row.names = 1))

head(insitu_matrix)
rownames(insitu_matrix[1:10,])
dim(insitu_matrix)

counts <- sobj_CPA@assays$RNA@counts %>%
  as.matrix() %>%
  t()
str(counts)
counts[25:30, 9:14]

data(ioprofiles)

# TODO: load safeTME cell profile matrix
# 2. generate cell profile matrix for tumor cells in the CPA

str(ioprofiles)
ioprofiles[1:5, 1:10]

negmean <- sobj_CPA@assays$negprobes@counts %>%
  as.matrix() %>%
  t() %>%
  Matrix::rowMeans()
head(negmean)

# sup <- insitutypeML(x = counts,
#                     neg = negmean,
#                     cohort = NULL,
#                     reference_profiles = insitu_matrix)   
# 
# sup$clus[1:10]
# 
# # nice plots
# p <- heatmap(sweep(sup$profiles, 1, pmax(apply(sup$profiles, 1, max), .2), "/"), scale="none",
#              main = "InSituType Clustering"
# )
# show(p)
# ggsave("Azimuth_heatmap_CPA.png", plot = p)
# 
# iocolors # list of hex codes
# cols <- iocolors[seq_along(unique(sup$clust))]
# names(cols) <- unique(sup$clust)
# 
# # object?
# plot(sobj_CPA@reductions$UMAP@cell.embeddings, pch=16, cex=0.75,
#      asp=1, cex.main=0.75, 
#      main = "cells in UMAP space",
#      col=cols[sup$clust],
#      xlab="", ylab="", xast="n", yast="n")
# legend("bottomleft", pch=16, col=cols, legend=names(cols), cex=0.7)
# 
# ggsave("Azimuth_CellsIn_umap_CPA.png")
# 
# 
# fp <- flightpath_plot(flightpath_result = NULL,
#                       insitutype_result = sup,
#                       col = cols[sup$clust])
# print(fp)
# ggsave("ImageDimPlot_AzCPA_FlightPath_CPA.png", plot = fp)
# 
# length(sup$clust)
# sobj_CPA <- AddMetaData(sobj_CPA, metadata = data.frame(CellType_AzCPA = sup$clust))
# 
# DefaultAssay(sobj_CPA) <- "RNA"
# Idents(sobj_CPA) <- 'CellType_AzCPA'
# #DefaultAssay(sobj_CPA) <- "RNA_normalized_cacc3861.a7b3.47a9.90f1.795a3d4cdcad_1"
# p <- ImageDimPlot(sobj_CPA, axes = TRUE, cols = "glasbey", dark.background = FALSE) 
# show(p)
# ggsave("ImageDimPlot_AzCPA_celltypes_CPA.png", plot = p)

###############################################
# Split Patient sobj and create FOVs
##############################################

## TODO: sort out PSR-01 and PSR-02 FOVs
# Split further
sobj_Org <- subset(sobj123, subset = Run_Tissue_name == 'PSR-Organoids')
sobj1 <- subset(sobj123, subset = Run_Tissue_name == 'PSR-01')
sobj2 <- subset(sobj123, subset = Run_Tissue_name == 'PSR-02')

# PSR-01 
sobj1@meta.data$imagerow <- sobj1@meta.data$CenterY_global_px
sobj1@meta.data$imagecol <- sobj1@meta.data$CenterX_global_px
# Ensure that spatial coordinates are set correctly
coords_df <- sobj1@meta.data[sobj1@meta.data$Run_Tissue_name == "PSR-01", c("imagerow", "imagecol")]

fov_psr01 <- CreateFOV(coords = coords_df, type = "centroids", assay="RNA", key="PSR-01")
sobj1@images$image <- fov_psr01
DefaultAssay(sobj1@images$image) <- "RNA"

Layers(sobj1)
Assays(sobj1)
DefaultAssay(sobj1) <- "RNA"
Idents(sobj1) <- 'Sample.ID'
#DefaultAssay(sobj_CPA) <- "RNA_normalized_cacc3861.a7b3.47a9.90f1.795a3d4cdcad_1"
p <- ImageDimPlot(sobj1, axes = TRUE, cols = "glasbey") 
# No compatible spatial coordinates present means a FOV (type SpatialImage) is needed. 
show(p)

# PSR-02
sobj2@meta.data$imagerow <- sobj2@meta.data$CenterY_global_px
sobj2@meta.data$imagecol <- sobj2@meta.data$CenterX_global_px
# Ensure that spatial coordinates are set correctly
coords_df <- sobj2@meta.data[sobj2@meta.data$Run_Tissue_name == "PSR-02", c("imagerow", "imagecol")]
coords_df[is.na(coords_df)] <- 0


fov_psr02 <- CreateFOV(coords = coords_df, type = "centroids", assay="RNA", key="PSR-02")
sobj2@images$image <- fov_psr02
DefaultAssay(sobj2@images$image) <- "RNA"

Layers(sobj2)
Assays(sobj2)
DefaultAssay(sobj2) <- "RNA"
Idents(sobj2) <- 'Sample.ID'
#DefaultAssay(sobj_CPA) <- "RNA_normalized_cacc3861.a7b3.47a9.90f1.795a3d4cdcad_1"
p <- ImageDimPlot(sobj2, axes = TRUE, cols = "glasbey") 
# No compatible spatial coordinates present means a FOV (type SpatialImage) is needed. 
show(p)

# many missing coords?
df_metadata <- sobj123@meta.data
coords_df <- df_metadata[df_metadata$Sample.ID == "ORG_2_3", c("Sample.ID", "CenterY_global_px", "CenterX_global_px")]
head(coords_df) # NA
sum(is.na(coords_df$CenterY_global_px)) # lots
sum(!is.na(coords_df$CenterY_global_px))

dim(coords_df)
fov_6_5 <- CreateFOV(coords = coords_df, type = "centroids", assay="RNA", key="PSR-02")

# or provide molecules. 

sobj2@images$image <- fov_6_5
sobj123@images



#####################################
# Cell Typing Plots
#####################################

sobj1 <- subset(sobj, subset = Run_Tissue_name == "PSR-01")

#Idents(sobj1) <- sobj1$Sample.ID
Idents(sobj1) <- 'Sample.ID'
p <- DimPlot(sobj1)
show(p)
ggsave("DimPlot_Azimuth_celltypes_PSR01.png", plot = p) 
p <- ImageDimPlot(sobj1, axes = TRUE, cols = "glasbey", dark.background = FALSE) 
show(p)
ggsave("ImageDimPlot_Azimuth_celltypes_PSR01.png", plot = p)

# PSR02 needs coords/FOV also

sobj2 <- subset(sobj, subset = Run_Tissue_name == "PSR-02")
#Idents(sobj1) <- sobj1$Sample.ID
Idents(sobj2) <- 'CellType'
p <- DimPlot(sobj2)
show(p)
#ggsave("DimPlot_Slide1.png", plot = p) 
p <- ImageDimPlot(sobj2, axes = TRUE, cols = "glasbey", dark.background = FALSE)
show(p)
ggsave("ImageDimPlot_Azimuth_celltypes_PSR02.png", plot = p)





sample1 <- subset(sobj1, subset = sample == "3_1")


#########################################
# Gene Aliases (move to utils)
########################################

genes_not_in_cosMx <- function(profile_genes, cosMx_genes){
  # Given 2 lists of genes, return the profile_genes not found in the cosMx list
  # for example, when checking a Cell Profile Matrix, we would like < 10% of genes to be mismatched. 
  missing <- list(setdiff(profile_genes, cosMx_genes))
  return(missing)
}

# look_for_aliases <- function(gene_list){
#   # given a list of genes, look up any aliases and return them in a nested list. 
#   
#   for gene in gene_list:
#     gene_aliases <- list()
#     gene_aliases[gene] <- list()
#     gene_aliases[gene] <- get_aliases(gene)
#     
#     
#     return(gene_aliases)
# }

# Function to look up gene aliases
lookup_gene_aliases <- function(gene, alias_list) {
  if (gene %in% names(alias_list)) {
    return(alias_list[[gene]])
  } else {
    return(NULL)
  }
}


library(biomaRt)

# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve HGNC symbols and their aliases
gene_aliases <- getBM(
  attributes = c("hgnc_symbol", "external_synonym"),
  mart = ensembl
)
class(gene_aliases)
print("BioMaRT")
gene_aliases[gene_aliases$hgnc_symbol == "NKX2-1", "external_synonym"]

# Install and load the hgnc package if not already installed
# if (!requireNamespace("hgnc", quietly = TRUE)) {
#   install.packages("hgnc")
# }

# library(hgnc)
# 
# # Query HGNC database for aliases
# hgnc_dataset <- import_hgnc_dataset(hgnc::latest_archive_url())
# alias_list <- hgnc_dataset[hgnc_dataset$symbol == "NKX2-1", "alias_symbol"] %>% pull(alias_symbol)
# alias_values <- alias_list[[1]]
# print("HGNC")
# alias_values


##########################################
#. Plots
#########################################
# cowplot, patchwork 

# how to make a grid of plots
library(patchwork)
# Create a list of ggplot objects
plots <- lapply(1:6, function(i) {
  ggplot(data.frame(x = 1:10, y = rnorm(10, i, 1)), aes(x = x, y = y)) +
    geom_point() +
    ggtitle(paste("Plot", i))
})

# Combine the plots
plot_grid <- plots[[1]] + plots[[2]] + plots[[3]] +
  plots[[4]] + plots[[5]] + plots[[6]] +
  plot_layout(ncol = 3, nrow = 2)

# Print the plot grid
plot_grid



####################################
# Giotto conversion, experimental; not used
####################################
# compare seuratToGiotto to createGiottoObject w raw files
# this will use > 300GB RAM, take > 1 hr and then crash R unless you limit the scope to a single slide. lung or PSR01 only is fast
# sobj1 <- subset(sobj, subset = Run_Tissue_name == "PSR-01")
df_newmetadata <- sobj1@metadata
cols_to_keep <- colnames(sobj1@meta.data[])
# make sure cell_ID dn exist in metadata
#[39] "cell_ID"    
sobj1[['cell_ID']] <- NULL
colnames(sobj1@meta.data)

# Caution: this function creates a column named "cell_ID" so make sure you don't have one already.
g1 <- seuratToGiottoV5(sobj = sobj1, 
                       spatial_assay = "SCT",  # or "Nanostring" 
                       #dim_reduction = c("pca", "umap"),
                       subcellular_assay = "Nanostring",
                       #sp_network = "knn",
                       #nn_network = "knn",
                       verbose=TRUE
)


# suspicious output imho:
# Selecting col "molecule" as feat_ID column
# Selecting cols "x" and "y" as x and y respectively
# 
# no external python path was provided, but a giotto python environment was found
# and will be used
# Consider to install these (optional) packages to run all possible Giotto
# commands for spatial analyses: trendsceek multinet RTriangle
# Giotto does not automatically install all these packages as they are not
# absolutely required and this reduces the number of dependencies
# Warning message:
#   Item 1 has 950 rows but longest item has 19001; recycled with remainder.

# Warning message:
#Item 1 has 950 rows but longest item has 19001; recycled with remainder. 


showGiottoFeatInfo(g1)
df_feat <- fDataDT(g1) # appears to have 950 features x the metadata. (?)
getFeatureMetadata(g1)
# it matched genes with cell meta; lets drop nearly all 

# Removing feature metadata this way caused problems. metadata referenced does not exist
# dim(df_feat)
# cols_to_drop <- colnames(df_feat)[!colnames(df_feat) %in% c("gene", "gene_ID", "feat_type", "feat_ID")]
# cols_to_drop[1]
# 
# for (i in 1:length(cols_to_drop)) {
#   g1@feat_metadata[cols_to_drop[i]] <- NULL
# }
# g1@feat_metadata['source'] <- 'seuratToGiottoV5'
# TODO: expand feat meta with GO terms, etc.
# also why are there repeats


showGiottoSpatialInfo(g1)
class(g1@cell_metadata$cell$rna) # cellMetaObj. rna vs negprobe 
#df_meta <- g1@cell_metadata$cell$rna@metaDT
df_meta <- pDataDT(g1)
g1@cell_metadata$cell$rna@spat_unit # cell
g1@cell_metadata$cell$rna@col_desc
g1@cell_metadata$cell$rna@feat_type # rna

g1@cell_metadata$cell$rna
colnames(df_meta)
#rownames(g1@cell_metadata$cell$rna@metaDT) <- g1@cell_metadata$cell$rna@metaDT$cell_ID

showGiottoCellMetadata(g1)

g1@spatial_locs # class spatLocsObj 
# sdimx  sdimy cell_ID
# <num>  <num>  <char>
#   1: 60649 129828 c_3_1_1

g1@spatial_locs[1] # has slots name, coordinates, spat_unit 
g1@spatial_locs$cell$raw@coordinates
g1@spatial_locs$cell$raw@spat_unit # "cell"

g1@spatial_info # class SpatVector, geometry: points, dimensions: 19001, 0 (geometries, attributes)
# SpatVector can represent points, lines, or polygons 
g1@spatial_info$cell # giottoPolygon class 
# giottoPolygon class has slots name, spatVector = terra spatVector to store polygon shapes, spatVectorCentroids, overlaps, 
g1@spatial_info$cell@spat_unit # "cell"

# subset 
df_meta = pDataDT(g1)
colnames(df_meta)
# weird, 'cell_ID' in metadata twice, but matches cell and cell_id
#g1 <- removeCellAnnotation(g1, columns=c("cell_ID")) # removes the first instances; leaves the second
cell_IDs_to_keep = df_meta[fov== 34]$cell_id
g_lymph_34 = subsetGiotto(g1, cell_ids = cell_IDs_to_keep)

cell_IDs_to_keep = df_meta[fov== 27]$cell_id
g_lung_27 = subsetGiotto(g1, cell_ids = cell_IDs_to_keep)

# but then.. 
#In FUN(X[[i]], ...) : spat_unit: cell spatloc name: raw
#cell IDs in spatial locations are missing from spatial polygon info

# 
spatPlot2D(g1)
spatPlot2D(g_lung_27)

# what are our cell types with counts
table(df_meta$final_cell_type)




