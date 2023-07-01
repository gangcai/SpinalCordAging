library(Seurat)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 20000 * 1024^2)
####### get all rds files ######
files.list <- list.files("../../../analysis/","rds")
results_all <- data.frame()
k <- 0
for(filename in files.list){
	k <- k + 1
	#### loading the data ####
	cellType <- sub("_SeuratV5.rds","",filename)
	print(cellType)
	obj <- readRDS(paste0("../../../analysis/",filename))
	DefaultAssay(obj) <- "RNA"

	######## cell-level wilcoxon two group DEG comparison ########
	### settings ####
	min_cell <- 3 # minimal number of cells for each group

	#### add group information #####
	Idents(obj) <- "groups"
	clusters <- unique(obj@meta.data$seurat_clusters)
	groups <- c("Old","Young") # change this according to your group information

	#### cluster-level wilcoxon test ######
	obj <- JoinLayers(obj)
	results <- ""
	cells <- names(obj$groups)
	i <- 0
	for(cid in clusters){
	   print(cid)
	   group1 <- groups[1]
	   cell.g1 <- cells[obj$groups==group1 & obj$seurat_clusters == cid]
	   group2 <- groups[2]
	   cell.g2 <- cells[obj$groups == group2 & obj$seurat_clusters == cid]
	   if(length(cell.g1) > min_cell & length(cell.g2) > min_cell){
		   markers <- FindMarkers(obj,ident.1=cell.g1, # ident.1 vs ident.2
					  ident.2=cell.g2,
					  slot="data",
					  verbose=F,
					  assay="RNA",
					  logfc.threshold=0,
					  only.pos=F,
					  test.use="wilcox",
					  pseudocount.use=0.1)
		   cluster.info <- rep(cid,nrow(markers))
		   cmp_sample <- paste0(group1,"_vs_",group2)
		   cmp_info <- rep(cmp_sample,nrow(markers))
		   hgr <- cbind(cluster.info,markers)
		   hgr <- cbind(cmp_info,hgr)
		   genes <- rownames(hgr)
		   hgr <- cbind(genes,hgr)
		   i <- i+1
		   if(i == 1){
		     results <- hgr
		   }else{
		     results <- rbind(results,hgr)
		   }
	   }
	}

	results$cellType <- cellType
	if(k == 1){
		results_all <- results
	}else{
		results_all <- rbind(results_all,results)
	}
}
##### generate the output file #######
write.table(results_all,
	    file="SC_subpopulations_DEGs_WilCox.tsv",
	    sep="\t",
	    quote=F,
	    row.names=F,
	    col.names=T)
