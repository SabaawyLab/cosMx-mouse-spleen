# cosMx-mouse-spleen
2024 cosMx mouse spleen TMA


Best code is the Quarto file, mouse_spleen_RNA.qmd which includes mouse_seurat_utils.R for helper functions. 

Purpose of script: load the Seurat data object produced by the AtoMx pipeline, 
generate plots, evaluate rna quality controls. summarize and display what has 
already been done in AtoMx 1

1. Seurat LoadNanostring() method to load Seurat from flatfiles 

2\. Load Sample metadata per cell:

2a. Subtract negative control probes
2b. Apply Quality Control filters. 
2c. Normalization method: Seurat LogNormalize
2d. FOV and Sample.ID metadata assignments. 
2e. Subset to slide labelled "Sabaawy new core 09/13/2024 5" and only spleen samples

3. Optimize PCA, Umap, assess clusters
4. Data Exploratory Analysis, FeaturePlots 
5. Integration, since Homozygous samples 
6. Apply Cell Typing to the Seurat object with clustering and manual annotation
with use of FindMarkers and DotPlot

7. plots and visualizations, zoomed in to representative regions
8. Export various metrics, e.g. percent of cells positive by marker

