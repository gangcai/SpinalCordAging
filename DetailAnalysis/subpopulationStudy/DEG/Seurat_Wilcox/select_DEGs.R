library(dplyr)
data <- read.table("SC_subpopulations_DEGs_WilCox_annotated.tsv",header=T,sep="\t")
padj_cutoff <- 0.05
log2fc_cutoff <- log2(1.5)
data_f <- data %>% filter((p_val_adj < padj_cutoff) & (abs(avg_log2FC) > log2fc_cutoff))

data_f <- data_f %>% mutate(deg_type=ifelse(avg_log2FC >0, "UP","DOWN"))

write.table(data_f,file="SC_subpopulations_DEGs_WilCox_annotated_sig.tsv",row.names=F,col.names=T,sep="\t",quote=F)
deg_stat <- data_f %>% group_by(deg_type,cellType,sub_cellType) %>%
		summarise(DEG_Num=n())

write.table(deg_stat,file="SC_subpopulations_DEGs_summary.tsv",row.names=F,col.names=T,sep="\t",quote=F)

