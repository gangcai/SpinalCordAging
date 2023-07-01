library(dplyr)
p_adj_cutoff <- 0.05
logfc_cutoff <- log2(1.5)

project_id <- "SCAge"
data <- read.table(paste0("../",project_id,"_DEGs_WilCox.tsv"),header=T,sep="\t")
data.f <- data %>% filter(p_val_adj < p_adj_cutoff & abs(avg_log2FC) > logfc_cutoff)
write.table(data.f,file=paste0(project_id,"_DEGs_WilCox_Sig.tsv"),col.names=T,row.names=F,sep="\t",quote=F)

#data.f$cluster.info <- factor(as.numeric(data.f$cluster.info))
sum.d <- data.f %>% 
    group_by(cluster.info) %>%
    summarise(no_rows = length(cluster.info))
colnames(sum.d) <- c("cluster","DEG_Num")
write.table(sum.d,file=paste0(project_id,"_DEGs_WilCox_Sig_Summary.tsv"),
	    col.names=T,row.names=F,sep="\t",quote=F)
