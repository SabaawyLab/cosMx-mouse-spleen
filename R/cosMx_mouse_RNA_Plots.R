# cosMx_mouse_RNA_Plots.R
## 
## ---------------------------
##
## Script name: cosMx_mouse_RNA_Plots.R
##
## Purpose of script: # Helper functions to work with cosMx Seurat objects 
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

# source("R/gene_annotations.R")") 

## ---------------------------


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

# #root_dir <- '/Volumes/T7Shield/PSR-GEN-057' # Experiment specific folder structure
# root_dir <- '/Users/annmstrange/Documents/cosMx/PSR-GEN-057'
# root_ref_dir <- '/Volumes/T7Shield' # reference dir for files not specific to this Experiment. 
# exp_name <- 'PSR-GEN-057'
# slide_names <- c( 'PSRCPA', 'PSROrganoids', 'PSR02', 'PSR01')
# # The cell ids in these are labeled c_1, c_2, c_3, c_4 respectively in the seurat object and metadata files.
# slide_names
# slide_name <- 'PSR01' # select one to work with (can be arbitrary)
# # run_names is not needed but could be used to navigate the directory structure
# run_names <- c('20240329_210634_S4', '20240329_210634_S3', '20240329_210634_S1', '20240329_210634_S2')
# run_name <- '20240329_210634_S2' # matches slide_name
# 
# # Caution: There is a cell naming convention in the flatfile vs the exported Seurat Object that differ. 
# # After seurat::LoadNanostring(), we will rename the keys to match the Seurat object, like c_<slide_num>_<fov>_<cell_num> e.g. c_1_1_300
# # cell_label_names is not explicitly needed, but merge will prepend these to the rownames if uniqueness is needed. 
# # Our use of rename_keys makes the use of these cell_label_names unnecessary on merge.
# # cell_label_names <- c('c_1', 'c_2', 'c_3', 'c_4')
# 
# results_dir <- paste0(root_dir, "/rna_results")
# results_sub_wholeassay <- "/WholeAssay"
# results_sub_patient <- "/ByPatient"
# results_sub_region <- "/ByPatientRegion"
# results_sub_celltype <- "/ByCellType"
# results_sub_panck <- "/ByPanCK"
# results_sub_annot <- "/ByAnnotatedRegion"
# results_sub_annot_celltype <- "/ByAnnotAndCellType"
# 
# rna_root_dir <- paste0(root_dir, '/rna')
# rna_dir <- paste0(root_dir, '/rna/PSR01/20240329_210634_S2/AnalysisResults/d02t7ym5l3') 
# # This text file lists all the Probes
# rna_lookup <- paste0(rna_root_dir, '/', slide_name, '/', run_name, '/plex-d02t7ym5l3.txt')
# 
# metadata_dir <- paste0(root_dir,'/Metadata')
# df_sample <- read.csv(paste0(metadata_dir, '/Sample_Metadata_DeIdent.csv'))
# df_fov_meta <- read.csv(paste0(metadata_dir, '/FOV_Sample_Ids_RNA.csv'))
# 
# # This nanostring-provided file lists Probe names, gene names, and some annotations
# markers_dir <- paste0(root_dir, '/BioMarkers')
# ns_gene_data <- 'LBL-11178-03-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx'
# ns_gene_fn <- paste0(markers_dir, '/', ns_gene_data)
# df_gene_data <- get_nanostring_gene_annotations(ns_gene_fn)
# 
# # drop ns_annotations
# #df_gene_data <- df_gene_data %>% dplyr::select(-ns_annotations)
# 
# df_special_genes <- read.csv(paste0(markers_dir, '/GenesOfInterest_RNA.csv'))
# colnames(df_special_genes) <- c('Display_Name', 'OfInterest_Reason', "Notes")
# head(df_special_genes)
# 
# genes_special <- list(df_special_genes$Display_Name[1:47])
# # breakdown
# genes_special_bypass <- df_special_genes %>%
#   filter(OfInterest_Reason %in% c("Bypass", "bypass", "PI3 Pathway", "MET pathway", "MAPK")) %>%
#   pull(Display_Name)
# genes_special_sq <- df_special_genes %>%
#   filter(OfInterest_Reason %in% c("plasticity", "Gentles Adeno or SQ", "SCC", "Gentles SCC", "SCC w KRT5")) %>%
#   pull(Display_Name)
# genes_special_sq
# 
# genes_special_sclc <- df_special_genes %>%
#   filter(OfInterest_Reason %in% c("SCLC", "SCLC progression")) %>%
#   pull(Display_Name)
# 
# genes_adeno <- df_special_genes %>%
#   filter(OfInterest_Reason %in% c("Adeno", "Gentles Adeno or SQ", "Gentles Adeno")) %>%
#   pull(Display_Name)

genes_special_ne <- c("S100B" )
genes_p53 <- c("TP43", "RB1")

# don't bother plotting 
genes_low_tumor_expr <- c("KRT5", "KRT13")

# left outer join with df_gene_data by RNA1K to gene 
df_gene_data2 <- left_join(df_gene_data, df_special_genes, by ='Display_Name')
df_gene_data2$OfInterest_Reason[is.na(df_gene_data2$OfInterest_Reason)] <- " "
head(df_gene_data2)
nrow(df_gene_data2)
df_gene_data <- df_gene_data2

file.exists(results_dir)
file.exists(rna_root_dir) 
file.exists(rna_dir)
file.exists(rna_lookup)
file.exists(metadata_dir)
file.exists(paste0(metadata_dir, '/Sample_Metadata_DeIdent.csv'))
file.exists(paste0(metadata_dir, '/FOV_Sample_Ids_RNA.csv'))

#################################
# Resume here
#################################
saved_rds_dir <- paste0(root_dir, "/saved_objs")
file.exists(saved_rds_dir)
sobj_load <- readRDS(paste0(saved_rds_dir, "/latest_sobj_16Sep24.rds"))
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
iocolors <- c("B-cell"= "darkblue", 
              "endothelial"=    "#996633" ,
              "fibroblast" = "#999999" ,
              "macrophage" =  "#006600" ,
              "mast" = "#FFFF00" ,
              "mDC" = "#00FF00"  ,
              "monocyte" = "#33CC00"  ,
              "neutrophil" =  "#9966CC"  ,
              "NK" =    "grey10" ,
              "pDC" = "#00FFFF" ,
              "plasmablast" =   "#3399CC"  ,
              "T CD4 memory" =  "#FF0000"  ,
              "T CD4 naive" =   "#CC0000",
              "T CD8 memory" =      "#FF9900" ,
              "T CD8 naive" =      "#FF6633",
              "Treg" =    "#FF66FF" )

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

saved_rds_dir
#saveRDS(sobj_spleen, paste0(saved_rds_dir,"/sobj_spleen.rds"))
sobj_spleen <- readRDS(file = paste0(saved_rds_dir, "/sobj_spleen.rds"))

# subsets to get this GSEA
# sobj_tumor <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") &
#                        CellType_Labeled.1 == "Tumor" &
#                        TimePoint %in% c("pre", "post"))
# sobj_tumor_pairs <- subset(sobj_tumor, subset = (PatientID %in% c(3,4,6) & Sample.ID != "ORG_4_1")
#                            | Run_Tissue_name == "PSR-CPA")
# 
# sobj_pt3_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 3)
# sobj_pt4_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 4)
# sobj_pt6_tumor <- subset(sobj_tumor_pairs, subset = PatientID == 6)




#####################################
# Bar Plot of cell counts by Sample
####################################

df_meta <- sobj_spleen@meta.data

table(df_meta$Sample.Nm)
table(df_meta$Genotype)
table(df_meta$Run_Tissue_name)

sobj1 <- subset(sobj_spleen, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj2 <- subset(sobj_spleen, subset = Run_Tissue_name == "Sabaawy new core 09/13/2024 6")

####################################
# compare different conditions
####################################
#sobj_spleen <- sobj_spleen_even 
save_dir <- paste0(results_dir, "/volcanos")
file.exists(save_dir)
if(!file.exists(save_dir)){
  dir.create(save_dir)
}

# Do we have too many SCTransforms? two is okay (one from each slide)
#sobj_spleen@assays$SCT@SCTModel.list 

# create various pairings
sobj_WTvHomo5 <- subset(sobj_spleen_even5, subset = Genotype %in% c("Spleen WT","Spleen Homozygous") & Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_WTvHetero5 <- subset(sobj_spleen_even5, subset = Genotype %in% c("Spleen WT","Spleen Heterozygous")& Run_Tissue_name == "Sabaawy new core 09/13/2024 5")
sobj_HeterovHomo5 <- subset(sobj_spleen_even5, subset = Genotype %in% c("Spleen Heterozygous", "Spleen Homozygous")& Run_Tissue_name == "Sabaawy new core 09/13/2024 5")

# redo SCTTransform completely 
# DefaultAssay(sobj_WTvHomo) <- "Nanostring"
# sobj_WTvHomo <- SCTransform(sobj_WTvHomo, assay="Nanostring", verbose = TRUE)

# to compare e.g. Sample.ID vs rest, add calculated column for Sample.ID vs rest
sobj_spleen_even5$Sample.30.vs.Rest <- ifelse(sobj_spleen_even5$Sample.ID == "30", "Sample.30", "Rest")
table(sobj_spleen_even5$Sample.30.vs.Rest)
#sobj_spleen$Sample.31.vs.Rest <- ifelse(sobj_spleen_even5$Sample.ID == "31", "Sample.31", "Rest")
#table(sobj_spleen_even5$Sample.31.vs.Rest)
sobj_spleen_even5$Sample.32.vs.Rest <- ifelse(sobj_spleen_even5$Sample.ID == "32", "Sample.32", "Rest")
table(sobj_spleen_even5$Sample.32.vs.Rest)
sobj_spleen_even5$Sample.34.vs.Rest <- ifelse(sobj_spleen_even5$Sample.ID == "34", "Sample.34", "Rest")
table(sobj_spleen_even5$Sample.34.vs.Rest)
#sobj_spleen$Sample.35.vs.Rest <- ifelse(sobj_spleen_even5$Sample.ID == "35", "Sample.35", "Rest")
#table(sobj_spleen_even5$Sample.35.vs.Rest)



table(sobj_WTvHomo5@meta.data$Genotype)

DefaultAssay(sobj_WTvHomo5) <- "RNA"
Idents(sobj_WTvHomo5) <- "Genotype"

#sobj_WTvHomo <- PrepSCTFindMarkers(sobj_WTvHomo, assay = "SCT", verbose = TRUE)

df_wt_v_homo <- get_DEGs(sobj_WTvHomo5, "Genotype", "Spleen Homozygous", df_gene_data, "RNA") 
print(paste("There are",nrow(df_wt_v_homo[df_wt_v_homo$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Genotype WT vs Homo 5"))
p <- plot_volcano_degs(df_wt_v_homo, num_labels = 50, title="WT vs Homozygous slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10, use_repel = FALSE) #  , special_labels = genes_special, 


# with genes of interest only labeled,
# try lower_num_labels
# If adj_p_val = 0, the volcanol plot will just peg those genes at the top. 
#.  make wider?


df_wt_v_hetero <- get_DEGs(sobj_WTvHetero5, "Genotype", "Spleen Heterozygous", df_gene_data) 
print(paste("There are",nrow(df_wt_v_hetero[df_wt_v_hetero$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Genotype WT vs Hetero 5"))
p <- plot_volcano_degs(df_wt_v_hetero, num_labels = 50, title="WT vs Heterozygous slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10, # special_labels = genes_special,
                       use_repel = FALSE)

df_hetero_v_homo <- get_DEGs(sobj_HeterovHomo5, "Genotype", "Spleen Homozygous", df_gene_data) 
print(paste("There are",nrow(df_hetero_v_homo[df_hetero_v_homo$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Genotype Hetero vs Homo 5"))
p <- plot_volcano_degs(df_hetero_v_homo, num_labels = 50, title="Heterozygous v Homozygous slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10, # special_labels = genes_special,
                       use_repel = FALSE)

# samples vs rest

df_30_v_rest <- get_DEGs(sobj_spleen_even5, "Sample.30.vs.Rest", "Rest", df_gene_data) 
print(paste("There are",nrow(df_30_v_rest[df_30_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 30 vs Rest"))
# rename 'Significant' to 'HighFoldChangeLowPval'
p <- plot_volcano_degs(df_30_v_rest, num_labels = 50, title="Sample 30 vs the Rest_slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10, # special_labels = genes_special,
                       use_repel = FALSE)

# df_31_v_rest <- get_DEGs(sobj_spleen5, "Sample.31.vs.Rest", "Rest", df_gene_data) 
# print(paste("There are",nrow(df_31_v_rest[df_31_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 31 vs Rest"))
# # rename 'Significant' to 'HighFoldChangeLowPval'
# p <- plot_volcano_degs(df_31_v_rest, num_labels = 100, title="Sample 31 vs the Rest", save_dir, 
#                        plot_height=10, plot_width=10)

df_32_v_rest <- get_DEGs(sobj_spleen_even5, "Sample.32.vs.Rest", "Rest", df_gene_data) 
print(paste("There are",nrow(df_32_v_rest[df_32_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 32 vs Rest"))
# rename 'Significant' to 'HighFoldChangeLowPval'
p <- plot_volcano_degs(df_32_v_rest, num_labels = 50, title="Sample 32 vs the Rest slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10, # special_labels = genes_special,
                       use_repel = FALSE)

# df_33_v_rest <- get_DEGs(sobj_spleen5, "Sample.33.vs.Rest", "Rest", df_gene_data) 
# print(paste("There are",nrow(df_33_v_rest[df_33_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 33 vs Rest"))
# # rename 'Significant' to 'HighFoldChangeLowPval'
# p <- plot_volcano_degs(df_33_v_rest, num_labels = 100, title="Sample 33 vs the Rest", save_dir, 
#                        plot_height=10, plot_width=10)

df_34_v_rest <- get_DEGs(sobj_spleen_even5, "Sample.34.vs.Rest", "Rest", df_gene_data) 
print(paste("There are",nrow(df_34_v_rest[df_34_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 34 vs Rest"))
# rename 'Significant' to 'HighFoldChangeLowPval'
p <- plot_volcano_degs(df_34_v_rest, num_labels = 50, title="Sample 34 vs the Rest slide5_Top50nolines", save_dir, 
                       plot_height=10, plot_width=10,  #special_labels = genes_special,
                       use_repel = FALSE)

# df_35_v_rest <- get_DEGs(sobj_spleen5, "Sample.35.vs.Rest", "Rest", df_gene_data) 
# print(paste("There are",nrow(df_35_v_rest[df_35_v_rest$HighFoldChangeLowPval == "Y",]), "Enr DEGs for Sample 35 vs Rest"))
# # rename 'Significant' to 'HighFoldChangeLowPval'
# p <- plot_volcano_degs(df_35_v_rest, num_labels = 100, title="Sample 35 vs the Rest", save_dir, 
#                        plot_height=10, plot_width=10)
# 

# write_csv(df_path_3, paste0(results_dir, results_sub_patient, results_sub_annot, "/DEGs_PathAnnotRegions_CUTO59.csv"))

# More volcanos, orig.
p <- plot_volcano_degs(df_wt_v_homo, num_labels = 100, title="WT vs Homozygous slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=20, use_repel = FALSE)
p <- plot_volcano_degs(df_wt_v_hetero, num_labels = 100, title="WT vs Heterozygous slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=20, # special_labels = genes_special,
                       use_repel = FALSE)
p <- plot_volcano_degs(df_hetero_v_homo, num_labels = 100, title="Heterozygous v Homozygous slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=20, # special_labels = genes_special,
                       use_repel = FALSE)

p <- plot_volcano_degs(df_30_v_rest, num_labels = 100, title="Sample 30 vs the Rest_slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=20, # special_labels = genes_special,
                       use_repel = FALSE)
p <- plot_volcano_degs(df_32_v_rest, num_labels = 100, title="Sample 32 vs the Rest slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=10, # special_labels = genes_special,
                       use_repel = FALSE)
p <- plot_volcano_degs(df_34_v_rest, num_labels = 100, title="Sample 34 vs the Rest slide5_Top100nolines", save_dir, 
                       plot_height=10, plot_width=10,  #special_labels = genes_special,
                       use_repel = FALSE)


# collect dfs 
df_sample2 <- df_sample %>%
  filter(Organ == "Spleen" & Sample.ID >= 30)

head(df_sample2)

wb <- createWorkbook()

# List of dataframes, slide 5 only
dfs <- list(Summary = df_sample2, 
            WT_v_Homoz = df_wt_v_homo, 
            WT_v_Hetero = df_wt_v_hetero,
            Hetero_vs_Homoz = df_hetero_v_homo,
            Sample30_v_Rest = df_30_v_rest,
            #Sample31_v_Rest = df_31_v_rest,
            Sample32_v_Rest = df_32_v_rest,
            #Sample33_v_Rest = df_33_v_rest,
            Sample34_v_Rest = df_34_v_rest )
            #Sample35_v_Rest = df_35_v_rest )


# Write each dataframe to a separate sheet
for (name in names(dfs)) {
  addWorksheet(wb, name)                 # Add a new sheet named after the dataframe
  writeData(wb, sheet = name, dfs[[name]])  # Write the dataframe to the sheet
}

# Save the workbook
saveWorkbook(wb, file = paste0(save_dir, "/Spleen DE Genes_Slide5_v2.xlsx"), overwrite = TRUE)


########################
# Volcano
###########



# Volcano plot

# 
# top_genes <- df_path_groups %>%
#   top_n(30, wt = abs(avg_log2FC))
# 
# ggplot(df_path_groups, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significant)) +
#   geom_point(alpha = 0.6) +
#   scale_color_manual(values = c("Y" = "red", "N" = "gray")) +
#   theme_minimal() +
#   labs(title = "Volcano Plot Pre v Post Tmt",
#        x = "Log2 Fold Change",
#        y = "-Log10 Adjusted P-Value") +
#   theme(legend.title = element_blank())+
#   # Add labels for the top 10 highest fold change genes
#   ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
#                            size = 3, 
#                            box.padding = 0.3, 
#                            max.overlaps = Inf)


# looser signif
df_celltype_tumor_4_ad$HighFoldChangeLowPval <- ifelse(df_celltype_tumor_4_ad < 0.05 & abs(avg_log2FC) > 1.0, "Y", "N")
table(df_celltype_tumor_4_ad$HighFoldChangeLowPval)

#save_dir = paste0(results_dir, results_sub_region, results_sub_annot_celltype)
p <- plot_volcano_degs(df_celltype_tumor_4_ad, num_labels = 20, title="Pre v Post Tumor cells in Tumor Region Pt4", save_dir) 
p
ggsave(paste0(results_dir, "/Volcano_Pt4_Adeno.png"), width = 4, height = 4, dpi=300)

p <- plot_volcano_degs(df_celltype_tumor_6_ad, num_labels = 100, title="Pre v Post Tumor cells in Tumor Region Pt6 Adeno", save_dir) 
ggsave(paste0(results_dir, "/Volcano_Pt6_Adeno.png"), width = 5, height = 5, dpi=300)


#df_celltype_tumor_6
#save_dir <- paste0(results_dir, results_sub_patient, results_sub_annot_celltype)
p <- plot_volcano_degs(df_celltype_tumor_6, num_labels = 100, title="Pre v Post Tumor cells in Tumor Region Pt6", save_dir) 
ggsave(paste0(results_dir, "/Volcano_Pt6.png"), width = 5, height = 5, dpi=300)

#save_dir <- paste0(results_dir, results_sub_region, results_sub_annot_celltype)
p <- plot_volcano_degs(df_celltype_tumor_4_ad_v_sclc, num_labels = 80, title="Adeno v SCLC Tumor cells in annot Region Pt4", save_dir)

#save_dir <- paste0(results_dir, results_sub_patient, results_sub_annot)
p <- plot_volcano_degs(df_path_4, num_labels = 20, title="Pre v Post Tmt Pathology Annotated Patient4", save_dir)
p

p <- plot_volcano_degs(df_path_6, num_labels = 50, title="Pre v Post Tmt Pathology Annotated Patient6", save_dir)

#save_dir <- paste0(results_dir, results_sub_patient, results_sub_annot_celltype)
p1 <- plot_volcano_degs(df_ct_annot_groups6, num_labels = 20, title="Tumor Celltype in Tumor Region Pt6", save_dir)
p2 <- plot_volcano_degs(df_panck_6, num_labels = 20, title="Tumor PanCK in Patient6", save_dir)


Idents(sobj) <- "path_annot_name.1"
ImageDimPlot(sobj, fov="PSR01", molecules = c("IGHA1",  "MZB1"), alpha=0.5, dark.background = FALSE) + 
  labs(title = paste("Slide1, (clockwise from top left): pt2 (lung), pt3 (liver cnb), pt3 (LN), pt1 (lung, tiny)"))
ggsave(paste0(results_dir, "/PSR01_IGHA1_MZB1.png"), width = 12, height = 10, dpi=300)

ImageDimPlot(sobj, fov="PSR02", molecules = c("IGHA1",  "MZB1"), alpha=0.5, dark.background = FALSE) + 
  labs(title = paste("Slide2, (clockwise from top left): pt5 (pleura), pt7 (T2), pt4 (liver), pt6 (LN), pt6 (brain), pt4 (liver)"))
ggsave(paste0(results_dir, "/PSR02_IGHA1_MZB1.png"), width = 12, height = 10, dpi=300)


#####################
# Pseudo bulk
##################  


sobj_panck_tumor <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") &
                             PanCK.PT == 1 &
                             TimePoint %in% c("pre", "post"))

sobj_cell_type_tumor <- subset(sobj, subset = Run_Tissue_name %in% c("PSR-01", "PSR-02") &
                                 CellType_Labeled.1 == "Tumor" &
                                 TimePoint %in% c("pre", "post"))

sobj_tumor <- sobj_panck_tumor

Idents(sobj_tumor) <- "TimePoint"

genes_special <- unlist(genes_special)
genes_w_expr <- c("CD9","EPCAM","KRT7","MET","GPNMB","S100A8","AGR2","AREG","CEACAM6")

p <- VlnPlot2(sobj_tumor, features = genes_w_expr, violin=F, pt=FALSE, group.by= "Sample.Label3", ncol=1)
p
ggsave(paste0(results_dir, "/genes_wexpr_vlnplot_tumor.png"), width=10, height=20,  dpi=300)

table(sobj@meta.data$PanCK.PT)
# sobj_tumor <- subset(sobj, subset = (PanCK.PT == 1 & Run_Tissue_name %in% c("PSR-01", "PSR-02"))
#                            & TimePoint %in% c("pre", "post") )
# 
# Extract the counts matrix and metadata
counts_matrix <- as.matrix(GetAssayData(sobj_tumor, layer = "counts"))  # Use "data" slot for normalized data if needed
sample_info <- sobj_tumor@meta.data$Sample.ID

# Create pseudo-bulk count matrix by aggregating tumor cells by sample (summing counts)
pseudo_bulk <- as.data.frame(counts_matrix) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "cell", value = "count", -gene) %>%
  left_join(data.frame(cell = colnames(counts_matrix), sample = sample_info), by = "cell") %>%
  group_by(gene, sample) %>%
  summarise(pseudo_bulk_counts = sum(count)) %>%
  spread(key = "sample", value = "pseudo_bulk_counts", fill = 0)

pseudo_bulk <- as.data.frame(pseudo_bulk)

rownames(pseudo_bulk) <- pseudo_bulk$gene
head(pseudo_bulk)
# Now 'pseudo_bulk' is a matrix where rows are genes and columns are samples, with aggregated counts for tumor cells.

# Load DESeq2
#library(DESeq2)

# Prepare DESeq2 dataset
# Assuming you have a metadata dataframe that includes your sample conditions (e.g., "Condition" column with "Pre-treatment" and "Post-treatment")
# sample_metadata <- data.frame(
#   sample = colnames(pseudo_bulk[,-1]),  # Sample IDs
#   condition = c(rep("Pre-treatment", 4), rep("Post-treatment", 5))  # Example: adjust according to your data
# )
# sample_metadata

# Example of metadata
# metadata <- data.frame(
#   sampleID = c("sample1", "sample2", "sample3", "sample4"),
#   condition = c("pre", "pre", "post", "post"),
#   row.names = "sampleID"
# )

# then get aggregate counts for each sample 
df_bulk_samples <- unique(sobj_tumor@meta.data[,c("Sample.ID", "TimePoint")])
#df_bulk_samples <- unique(sobj_tumor@meta.data[,"TimePoint"])
colnames(df_bulk_samples) <- c("sample", "condition")
row.names(df_bulk_samples) = df_bulk_samples$sample
df_bulk_samples$condition <- factor(df_bulk_samples$condition)
#colnames(df_bulk_samples) <- "condition"
# resort by rownames
df_bulk_samples <- df_bulk_samples[order(rownames(df_bulk_samples)), ]

head(df_bulk_samples)

# check ncol(countData) == nrow(colData)
ncol(pseudo_bulk[,-1]) # 2 pre vs post
nrow(df_bulk_samples) # 9 sample ids

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_bulk[,-1],  # Exclude the gene column
  colData = df_bulk_samples,
  design = ~ condition  # Specify your experimental design
)

# Run DESeq2 for differential expression analysis
dds <- DESeq(dds)
results <- results(dds)

results$Gene <- rownames(results)
# move Gene to first column
results <- results[,c(ncol(results), 1:(ncol(results)-1))]
# replace NA with 0 for padj
results$padj[is.na(results$padj)] <- 1
results$padj[is.na(results$log2FoldChange)] <- 0

# View top differentially expressed genes
head(results)


# save results
write.csv(as.data.frame(results), file = paste0(results_dir, "/DESeq2_panck_results_16Sep24.csv"),
          row.names = TRUE)

# Volcano Plots

# calc Signif
results <- as.data.frame(results)
results$HighFoldChangeLowPval <- ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 0.6, "Y", "N")

top_genes <- as.data.frame(results) %>%
  dplyr::top_n(30, wt = abs(log2FoldChange))

ggplot(results, aes(x = log2FoldChange, y = -log10(padj) , color = HighFoldChangeLowPval)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Y" = "red", "N" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot Pre v Post Tmt",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank()) +
  # Add labels for the top 10 highest fold change genes
  ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                           size = 3, 
                           box.padding = 0.3, 
                           max.overlaps = Inf)


p <- plot_deseq_volcano_degs(results, num_labels = 50, title="Pre v Post Tmt CellType Tumor", results_dir) 
p
ggsave(paste0(results_dir, "/volcano_DESeq_pre_v_post_tmt_celltype_tumor.png"), p, dpi=300)




####################################
# GSEA
####################  



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




