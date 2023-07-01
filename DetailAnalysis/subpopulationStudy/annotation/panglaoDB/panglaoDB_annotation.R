library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
pdb <- read.table("PanglaoDB_markers_27_Mar_2020.tsv",header=T,sep="\t")
cluster_marker_df <- read.table("../../analysis/Microglia_merged_markers.tsv",header=T,sep="\t")

cluster_marker_df$geneU <- sapply(cluster_marker_df$gene, toupper)

celltypes <- unique(pdb$cell.type)
clusters <- unique(cluster_marker_df$cluster)

results <- data.frame()
k <- 0
for(cid in clusters){
  k <- k + 1
  marker_gene <- cluster_marker_df[cluster_marker_df$cluster == cid,]$geneU
  cname <- paste0("cluster",cid)
  
  item <- ""
  i <- 0
  for(celltype in celltypes){
    pdb_gene <- pdb[pdb$cell.type == celltype,]$official.gene.symbol
    total <- length(pdb_gene)
    overlapped_genes <- pdb_gene[pdb_gene %in% marker_gene]
    overlapped_num <- length(overlapped_genes)
    
    o_genes <- toString(overlapped_genes)
    ratio <- overlapped_num/total
    i <- i + 1
    info <- c(cname, celltype, ratio, overlapped_num, total,length(marker_gene),o_genes)
    if(i == 1){
      item <- info
    }else{
      item  <- rbind(item , info)
    }
  }
  colnames(item) <- c("cluster_id","celltype","pdb_marker_ratio","overlapped_num","pdb_total","cluster_marker_num","overlapped_genes")
  
  item[,"pdb_marker_ratio"] <- as.numeric(item[,"pdb_marker_ratio"])
  
  item <- item[order(item[,"pdb_marker_ratio"],decreasing = T),]
  if(k == 1){
    results <- item
  }else(
    results <- rbind(results,item)
  )
}

results <- as.data.frame(results)

results_top <- results %>% group_by(cluster_id) %>% top_n(10,pdb_marker_ratio)
results_top$pdb_marker_ratio <- as.numeric(results_top$pdb_marker_ratio)
results_top$overlapped_num <- as.numeric(results_top$overlapped_num)

plot.list <- list()
i <- 0
for(cid in unique(results_top$cluster_id)){
  i <- i + 1
  plot.list[[i]] <- ggplot(results_top[results_top$cluster_id == cid,], aes(x=pdb_marker_ratio,y=celltype, size=overlapped_num, col=pdb_marker_ratio)) + 
    geom_point() +
    ggtitle(cid)
}
p <- wrap_plots(plot.list,ncol=3)
pdf("panglaoDB_enrichment.pdf",width = 20,height=20)
print(p)
dev.off()

write.table(results, file="panglaoDB_enrichment_results.tsv",sep="\t",row.names=F,col.names=T,quote=F)
