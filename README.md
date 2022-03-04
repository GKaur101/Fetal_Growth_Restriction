Code for FGR/IUGR analyses

Seurat objects used in these scripts are located on the SCP:
https://singlecell.broadinstitute.org/single_cell/study/SCP1312/
In the file: IUGR_SeuratObjects.zip
- gcdata_02162018_ss2_NK_Seurat3.Rda
- gcdata_02072018_10x_TD_Seurat3.Rda
- gcdata_08072018_10x_NK_Seurat3.Rda

github_DE_genes_IUGR.Rmd
- this script generates DE genes between IUGR, CTL1, CTL2 (pairwise) using pseudobulk and mixed effects model approaches
- calls de_genes_mixed_effects_IUGR.R
- calls de_genes_pseudobulk_IUGR.R

github_topic_models_IUGR.Rmd
- this script generates a topic model for k = 16 and tol = 0.01, the final chosen parameters for our IUGR topic modeling analyses
- calls plot.umap.feature_seurat3.R 

