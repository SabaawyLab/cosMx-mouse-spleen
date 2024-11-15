# helper functions
# some might not work, but they are here for reference



# Basic function to convert human to mouse gene names using BiomaRt
convertHumanGeneList <- function(human_genes){
  # set host to "dec2021.archive.ensembl.org" to use the archived version of Ensembl
  # workaround for bug in latest host: https://support.bioconductor.org/p/9143914/, set host to prior version 105.
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "https://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://dec2021.archive.ensembl.org")
  
  lookup_list <- list()
  for (human_gene in human_genes){
    
    lookup_pair <- getLDS(attributes = c("hgnc_symbol"),
                          filters = "hgnc_symbol", values = human_gene, mart = human,
                          attributesL = c("mgi_symbol"), martL = mouse)
    
    if (dim(lookup_pair)[1] == 0){
      print(paste("No match found for", human_gene))
      lookup_list[human_gene] <- ""
    } else {
      lookup_list[human_gene] <- unique(lookup_pair[1, 2])
    }
    
    print(lookup_pair[,2])
    
    #lookup_list[[human_gene]] <- unique(genesV2[, 2])
  }
  
  #mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(mousex))
  return(lookup_list)
}



plot_waterfall_degs <- function(sobj, df_degs, num_labels = 50, title="", results_dir, 
                                plot_height=10, plot_width=10, special_labels=NULL, use_repel=TRUE) {
  # Expect a column names "HighFoldChangeLowPval" with values "Y" or "N"
  # top_genes <- df_degs %>%
  #   dplyr::filter(HighFoldChangeLowPval == "Y") 
  
  if (missing(special_labels)) {
    special_labels <- c()
    print ("No special labels, proceed with num_labels")
    top_genes <- df_degs %>%
      dplyr::top_n(num_labels, wt = abs(avg_log2FC)) %>%
      dplyr::select(1:10)
    print(head(top_genes))
  } else {
    top_genes <- df_degs %>%
      dplyr::filter(Gene %in% special_labels & HighFoldChangeLowPval == "Y") %>%
      dplyr::select(1:10)
    #top_genes <- as.data.frame(special_labels)
    df_degs$HighFoldChangeLowPval <- ifelse(df_degs$Gene %in% special_labels & 
                                              df_degs$HighFoldChangeLowPval == "Y", "Special", 
                                            df_degs$HighFoldChangeLowPval)
    print(head(top_genes))
  }
  
  # Cap very low p-values at a threshold (e.g.,tried 1e-300, now 1e-400) to avoid infinite values
  df_degs$p_val_adj <- pmax(df_degs$p_val_adj, 1e-300)
  
  p <- WaterfallPlot(sobj, features = top_genes)
  
  print(p)
  #ggplot2::ggsave(paste0(results_dir, "/WaterfallPlot_", sanitize_name(title), ".png"), 
  #                plot=p, 
  #                dpi = 300)
  return(p)
}

# fetch_gene_data <- function(sobj, genes_list){
#   # given a subset object and list of genes, fetch normalized data
#   # and metadata needed for typical boxplots
#   norm_data <- FetchData(sobj, gene_names, layer = "data")
#   df_meta <- sobj@meta.data[c('TimePoint','Sample.ID','Sample.Label', 'PatientID')]
#   gene_data <- merge(norm_data, df_meta, how="inner", by="row.names")
#   return (gene_data)
#   
# }
# 
# 

################################
# sub-fovs, Future use, alternate way to get a new fov to assign an sobj
################################

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

################################
# Inspect Sample 6 (brain)
################################


# sobj6 <- subset (sobj, subset = Sample.ID %in% c(6))
# df_meta6 <- sobj6@meta.data
# bbox_mtx <- get_bbox_of_sample(df_meta6, "Sabaawy new core 09/13/2024 5", 6)
# bbox_mtx
# cropped.coords <- Crop(sobj[["Sabaawynewcore091320245"]], 
#                        x = c(bbox_mtx[2,1], bbox_mtx[2,2]), y = c(bbox_mtx[1,1], bbox_mtx[1,2]), coords = "plot")
# sobj[["fov6"]] <- cropped.coords
# DefaultBoundary(sobj[["fov6"]]) <- "centroids"
# genes_regions6 <- c("Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd8b1", "Foxp3")
# colors_regions6 <-  c("red", "red", "red", "green", "green", "blue")
# names(colors_regions6) <- genes_regions6
# 
# p <- plot_molecules (sobj, "fov6", mol_list = genes_regions6, 
#                 colors_regions6, "Regions Sample 6", show_legend=TRUE, molsize=1)
# p
# ggsave(paste0(results_dir,"regions_sample6.png"), p, width=10, height=10)
# 
# # cell counts by fov
# count_mtx6 <- FetchData(sobj, fov = "fov6", vars=genes_regions6, layer = "counts")
# head(count_mtx6)
# 
# file_suffix <- "Brain6"
# gene_data1 <- fetch_mouse_gene_data(sobj6, gene_names6)
# head(gene_data1)
# p <- plot_mouse_multiple_gene_boxplots(gene_data1, gene_names)+
#   ggtitle("Sample 6")
# p
# 
# DefaultAssay(sobj6) <- "Nanostring"
# counts_data <- FetchData(sobj6, genes_regions6, layer = "counts") # s/b data
# df_meta6 <- sobj6@meta.data[c('Sample.ID','Sample.Label', 'Patient', 'fov')]
# head(df_meta6)
# gene_data <- merge(counts_data, df_meta6, how="inner", by="row.names")
# plot_mouse_multiple_gene_boxplots(gene_data, genes_regions6)+
#   ggtitle("Sample 6")
# 
# head(gene_data)

# I want to plot the number of rows with Cd3e > 0, Cd3d > 0, Cd3g > 0, or Cd8a > 0 
# for each fov
# "Cd8a", "Cd8b1"

# summary_df <- gene_data %>%
#   filter(Foxp3 > 0 ) %>%  # Filter rows where any of the conditions are true
#   group_by(fov) %>%  # Group by fov
#   summarise(FoxP3_pos_cells = n())  
# 
# summary_df


