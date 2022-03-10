library(Seurat)
library(dplyr)
library(patchwork)
source("ReadData.R")


UMAP <- function(){
  path <- paste(data_path, "/initial_data.rds", sep = "")
  if(file.exists(path)){
    dat <- readRDS(path)
  }else{
    data_path <- paste(getwd(),"/data/hg19",sep = "")
    dat <- ReadData(data_path)}
    # UMAP
  dat <- RunUMAP(dat, dims = 1:10, label = T)
  cluster_ids_path <- paste(data_path,"/cluster_ids.txt",sep = "")
  if(file.exists(cluster_ids_path)){
  new.cluster.ids <- c(read.table(cluster_ids_path,sep='\n')$V1)
  names(new.cluster.ids) <- levels(dat)
  dat <- RenameIdents(dat, new.cluster.ids)
  }
  return(dat)
}
data_path <- paste(getwd(),"/data/hg19",sep = "")
img_path <- paste(getwd(),"/img",sep = "")
marker_genes <- c(read.table(paste(data_path,"/marker_genes.txt",sep = ""),
                             sep = "\n")$V1)
dat <- UMAP()
p1 <- DimPlot(dat, reduction = 'umap', label = T, pt.size = 0.5) + NoLegend()
p2 <- FeaturePlot(dat, features = marker_genes)
svg(paste(img_path, "/UMAP.svg",sep= ""),width=10,height=20)
p1/p2
dev.off()
