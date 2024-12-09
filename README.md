# cosMx-mouse-spleen
2024 cosMx mouse spleen TMA


Best code is the Quarto file, mouse_spleen_RNA.qmd which includes mouse_seurat_utils.R for helper functions. 

What this Quarto file does:
1. read in mouse spleen flatfiles exported from AtoMx.
2. read in custom csv files with fov and sample metadata for mice. 
3. read in BioMarker of interest csv and Panel gene list
4. Integrationn to handle sample differences
5. standard Seurat steps
6. Cell Typing by clustering, and also InSituType (not used)
