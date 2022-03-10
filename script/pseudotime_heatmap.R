library(monocle)
library(dplyr)
library(patchwork)

data_path <- paste(getwd(),"/data/hg19",sep = "")
img_path <- paste(getwd(),"/img",sep = "")
rds_path <- paste(data_path,"/monocle_data.rds",sep = "")
HSMM <- readRDS(rds_path)
marker_genes <- c(read.table(paste(data_path,"/marker_genes.txt",sep = ""),
                             sep = "\n")$V1)
marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% marker_genes))
diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

svg(paste(img_path, "/pseudotime_heatmap.svg",sep= ""),width=10,height=10)
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)
dev.off()
