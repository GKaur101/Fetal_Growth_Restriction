## Written by Caroline Porter, November 18, 2017
## Colors function provided by Sam Riesenfield 
## This funciton takes a feature from the meta data of a Seurat object, or the un-scaled gene expression a 
## Seurat object, and plots it on UMAP with a blue to yellow to red gradient of colors. If more than one feature
## is supplied, the graphs will be sub-plotted. As it's written, meta data features and genes cannot be plotted
## simultaneously. 
## INPUTS
## obj - the Seurat object
## features - list of features to be plotted from meta data, or list of genes to plot
## feature.type - should be either "gene" or "meta" to indicate the type of feature
## ncols - number of columns desired for facet wrap
## pt.size - desired plot point size
## same.scale - TRUE if you would like your color bar scaled based on entire dataset rather than just the expression for plotted features 
## same.scale only makes sense when plotting genes
## OUTPUT
## a ggplot object is returned 

plot.umap.feature <- function(obj=countData, features="nGene", feature.type="meta", 
                              ncols=ceiling(sqrt(length(features))), pt.size=1, same.scale=FALSE, title="val", 
                              lower=NULL, upper=NULL, na.color="gray") {
        
# load required libraries 
library(reshape2)
library(plotly)
library(tidyr)

colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
"#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
"#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000")

# error out if user doesn't specify correct feature type 
if (feature.type!="meta" & feature.type!="gene"){stop("feature type must be 'meta' or 'gene'")}

# collect feature info
if (feature.type=="meta"){ 
        feature.info <- as.matrix(obj@meta.data[,features]) 
        if (length(features)==1){
                colnames(feature.info) <- features
                # feature.info <- t(feature.info)
        } 
        
        # build data frame of feature info and tsne coordinates 
        tmp.df <- data.frame(feature.info, obj@reductions$umap@cell.embeddings)
        plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
        
} else if (feature.type=="gene"){ 
        feature.info <- as.matrix(GetAssayData(object = obj, slot = "data")[features,]) 
        if (length(features)==1){
                colnames(feature.info) <- features
                feature.info <- t(feature.info)
        }   
        
        # build data frame of feature info and tsne coordinates 
        tmp.df <- data.frame(t(feature.info), obj@reductions$umap@cell.embeddings)
        plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
}



# set scales for color bar
if (!length(lower) & !length(upper)){
        lower = min(plot.df$val)
        upper = max(plot.df$val)
}



# plot the data on UMAP
p <- ggplot(plot.df, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=val), alpha=0.8, shape=19, size=pt.size) + 
        theme(aspect.ratio = 1) + scale_color_gradientn(colors=colors, limits=c(lower, upper), na.value=na.color)
p <- p + theme(aspect.ratio=1, text = element_text(size=10), axis.text=element_text(size=6), 
               strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))) + 
                labs(color = title) +
                facet_wrap( ~ name, ncol=ncols) + 
                theme(strip.text = element_text(size=10),
                      panel.background = element_blank(), 
                      plot.background = element_blank(),
                      axis.line = element_line(size=0.5))
return(p)
}


