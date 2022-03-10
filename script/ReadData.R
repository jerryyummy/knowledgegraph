library(Seurat)
library(dplyr)
library(patchwork)

ReadData <- function(data_path){
  path <- paste(data_path, "/initial_data.rds", sep = "")
  if(file.exists(path)){
    dat <- readRDS(path)
    return(dat)
  }else{
    dat.data <- Read10X(data_path)
    dat <- CreateSeuratObject(counts = dat.data, project = "Seurat", min.cells = 3, min.features = 200)
    # QC
    dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
    dat <- subset(dat, subset =  nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    # 标准化
    dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
    dat <- FindVariableFeatures(dat,selection.method = "vst", nfeatures = 2000)
    # 分类
    all.genes <- rownames(dat)
    dat <- ScaleData(dat, features = all.genes)
    dat <- RunPCA(dat, features = VariableFeatures(dat))
    dat <- JackStraw(dat, num.replicate = 100)
    dat <- ScoreJackStraw(dat, dims = 1:20)
    dat <- FindNeighbors(dat, dims = 1:10)
    dat <- FindClusters(dat, resolution = 0.5)
    saveRDS(dat,path)
    return(dat)
  }
}

