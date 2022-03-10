library(monocle)
library(dplyr)
library(patchwork)
source("ReadData.R")

data_path <- paste(getwd(),"/data/hg19",sep = "")
img_path <- paste(getwd(),"/img",sep = "")
dat <- ReadData(data_path)
data <- as(as.matrix(dat@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame',data=dat@meta.data)
fData <- data.frame(gene_short_name=row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame',data=fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
HSMM <- monocle_cds
# ¹éÒ»»¯
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >=10))
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 10,
                        reduction_method = 'tSNE', verbose = 'T')
HSMM <- clusterCells(HSMM, num_clusters = 5)
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Size_Factor + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 5)
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = '~percent.mt')
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
# ½µÎ¬
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

p1 <- plot_cell_trajectory(HSMM, color_by = 'seurat_clusters')
p2 <- plot_cell_trajectory(HSMM, color_by = 'State')
p3 <- plot_cell_trajectory(HSMM, color_by = 'Pseudotime')
svg(paste(img_path, "/pseudotime.svg",sep= ""),width=10,height=10)

p1/p2/p3
dev.off()


path <- paste(data_path, "/monocle_data.rds", sep = "")
saveRDS(HSMM,path)
