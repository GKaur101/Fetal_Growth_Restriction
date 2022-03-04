de_genes_mixed_effects_IUGR <- function(gcdata.cluster, cond1, cond2, filename) {

  # Subset cluster-specific Seurat object to just these two conditions 
  cond.data <- subset(gcdata.cluster,
                      cells=colnames(gcdata.cluster)[gcdata.cluster@meta.data$genotype==cond2 |
                                                       gcdata.cluster@meta.data$genotype==cond1])
  # Create data matrices
  data = GetAssayData(cond.data)
  data.2 <- data[,cond.data@meta.data$genotype==cond2, drop=F]
  data.1 <- data[,cond.data@meta.data$genotype==cond1, drop=F]
  
  # Find genes above min expression
  frac.1 <- rowSums(data.1 > 0)/dim(data.1)[2]
  frac.2 <- rowSums(data.2 >0)/dim(data.2)[2]
  genes.frac <- union(names(frac.1)[frac.1>=.1], names(frac.2)[frac.2>=.1])
  
  # Make a list of genes to test (genes expressed in >10% of cells in at least one condition)
  genes.test <- genes.frac[1:5]
  print(length(genes.test))
  
  # mixed-effects Poisson regression model
  library(lme4)
  f1= 'gene~condition+offset(log(scalenUMI))+(1|Sample)'
  results=NULL
  df=GetAssayData(cond.data, slot="counts") # collect counts 
  df=df[genes.test,] # subset to just genes being tested 
  
  # collect condition and covariate information (sample + UMI)
  meta=data.frame(cellnames=colnames(cond.data),condition=cond.data@meta.data$genotype,
                  Sample=cond.data@meta.data$mouse, nUMI=cond.data@meta.data$nCount_RNA)
  
  # loop through all genes being tested, run DE test, and store result  
  for (gene in c(1:dim(df)[1])){
    print(gene)
    dat=data.frame(gene=df[gene,],meta)
    dat$condition=factor(dat$condition,levels=c(cond2, cond1))
    mc=mean(dat$nUMI)
    dat=dat%>%mutate(scalenUMI=nUMI/mc)
    m <- glmer(f1, data =dat, family = poisson())
    coefs <- data.frame(coef(summary(m)))
    tmp=NULL
    for(j in c(1:2)){
      if (j==1){
        x=coefs[j,]
        names(x)=paste0(row.names(coefs)[j],"_",names(x))
        tmp=data.frame(celltype=i, gene=rownames(df)[gene],x) #
      } else{
        x=coefs[j,]
        names(x)=paste0(row.names(coefs)[j],"_",names(x))
        tmp=data.frame(tmp,x)}
    }
    results=rbind(results,tmp)
  }
  
  # format output 
  rownames(results)=results$gene
  colnames(results)=gsub("X.","",colnames(results))
  names(results)[10]="pval"
  
  # calculate FDR 
  results$FDR <- p.adjust(results$pval, method="fdr")
  
  # save results 
  write.csv(results, filename, row.names=F)
  
}