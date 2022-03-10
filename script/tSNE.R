library(Seurat)
library(dplyr)
library(patchwork)
source("ReadData.R")
data_path <- paste(getwd(),"/data/hg19",sep = "")
img_path <- paste(getwd(),"/img",sep = "")
dat <- ReadData(data_path)

# tSNE
top10 <- head(VariableFeatures(dat), 10)
dat <- RunTSNE(dat, dims = 1:10)
p1 <- DimPlot(dat, reduction = 'tsne')
p2 <- FeaturePlot(dat, features = top10)
svg(paste(img_path, "/tSNE.svg",sep= ""),width=10,height=10)
p1/p2
dev.off()
