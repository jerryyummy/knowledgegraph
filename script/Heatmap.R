library(Seurat)
library(dplyr)
library(patchwork)
source("UMAP.R")

data_path <- paste(getwd(),"/data/hg19",sep = "")
img_path <- paste(getwd(),"/img",sep = "")

dat <- UMAP()
marker_genes <- c(read.table(paste(data_path,"/marker_genes.txt",sep = ""),
                             sep = "\n")$V1)
# »æÍ¼
svg(paste(img_path, "/Heatmap.svg",sep= ""),width=10,height=10)
DoHeatmap(subset(dat, downsample = 100), features = marker_genes, size = 3)
dev.off()
