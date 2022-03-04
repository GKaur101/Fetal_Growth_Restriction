de_genes_pseudobulk_IUGR <- function(gcdata.cluster, cond1, cond2, save.path.sub, cluster, ext) {
  
    library(Matrix.utils)
    library(DESeq2)
    library(pheatmap)
    library(ggplot2)
    library(cowplot)
    
    # Subset to just the two conditions being tested 
    cond.data <- subset(gcdata.cluster,
                        cells=colnames(gcdata.cluster)[gcdata.cluster@meta.data$genotype==cond2 |
                                                                  gcdata.cluster@meta.data$genotype==cond1])
    # Split seurat object by mouse 
    gcdata.mouse <- SplitObject(cond.data, split.by="mouse")
    
    # Collect raw counts for each mouse 
    data.mouse <- lapply(names(gcdata.mouse), function(x) gcdata.mouse[[x]]@assays$RNA@counts)
    names(data.mouse) <- names(gcdata.mouse)
    
    # Pseudobulk - summing (rather than averaging) counts across all cells from a mouse according to 
    # https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
    sum.mouse <- lapply(names(data.mouse), function(x) rowSums(data.mouse[[x]]))
    names(sum.mouse) <- names(data.mouse)
    
    # collect meta data 
    meta <- data.frame(mouse=cond.data@meta.data$mouse, pheno=cond.data@meta.data$genotype)
    meta <- unique(meta)
    rownames(meta)<-meta$mouse
    meta$pheno <- factor(meta$pheno, levels=c(cond2, cond1))
    
    # Make counts matrix of all mice 
    x = do.call(cbind.data.frame, sum.mouse)
    
    # Only test genes that are expressed only in 10% of cells in one condition in this cluster 
    data = GetAssayData(cond.data)
    data.2 <- data[,cond.data@meta.data$genotype==cond2, drop=F]
    data.1 <- data[,cond.data@meta.data$genotype==cond1, drop=F]
    
    # Find genes above min expression
    frac.1 <- rowSums(data.1 > 0)/dim(data.1)[2]
    frac.2 <- rowSums(data.2 >0)/dim(data.2)[2]
    genes <- union(names(frac.1)[frac.1>=.1], names(frac.2)[frac.2>=.1])
    print(length(genes))
    
    # DESeq2
    ss2_DDS <- DESeqDataSetFromMatrix(countData = x[genes,], colData = meta, design = ~pheno)
    
    # PCA and heatmap 
    rld <- rlog(ss2_DDS, blind=TRUE)
    DESeq2::plotPCA(rld, intgroup = "pheno")
    ggsave(paste0(save.path.sub, "/DESeq2_PCA_cluster_", cluster, "_", cond2, "vs", cond1, ext, '.png'), height=4, width=4)
    
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    png()
    p = pheatmap(rld_cor, annotation = meta[, c("pheno"), drop=F])
    save_plot(paste0(save.path.sub, "/DESeq2_heatmap_cluster_", cluster, "_", cond2, "vs", cond1, ext, '.png'), p, base_width=4)
    dev.off()
    
    # Find DE genes 
    dds <- DESeq(ss2_DDS)
    plotDispEsts(dds)
    
    # Output results of Wald test for contrast for stim vs ctrl
    print(paste(levels(meta$pheno)[1], 'vs', levels(meta$pheno)[2]))
    
    contrast <- c("pheno", levels(meta$pheno)[1], levels(meta$pheno)[2])
    
    # resultsNames(dds)
    res <- results(dds, contrast = contrast)
    res <- lfcShrink(dds, contrast = contrast, res=res)
    
    # FDR using empirical null 
    # remove genes filtered out by independent filtering and the dispersion outliers
    res2 <- res[!is.na(res$padj), ]
    res2 <- res[!is.na(res$pvalue), ]
    
    # use z–scores returned by DESeq2 as input to fdrtool
    library(fdrtool)
    corrected <- fdrtool(res2$stat, statistic= "normal", plot = T)
    corrected_qval <- data.frame(qval=corrected$qval)
    rownames(corrected_qval) <- rownames(res2)
    corrected_qval$gene <- rownames(corrected_qval)
    corrected_qval_subset <- subset(corrected_qval, corrected_qval$qval<0.1)
    
    full <- cbind(res2, corrected_qval)
    colnames(full)[grep("qval", colnames(full))] <- "empirical_null_FDR"
    
    # small <- as.data.frame(subset(res, res$padj<0.1))
    write.csv(full, paste0(save.path.sub, "/DESeq2_cluster_", cluster, "_", cond2, "vs", cond1, ext, '.csv'))
}