library(ggplot2)
#library(tidyverse)
library(ggrepel)
library(dplyr)
p_adj_cutoff <- 0.05
logfc_cutoff <- log2(1.5)
top_num <- 5
data <- read.table("../SCAge_DEGs_WilCox.tsv",header=T,sep="\t")
data$label <- ifelse((data$p_val_adj < p_adj_cutoff & abs(data$avg_log2FC)> logfc_cutoff), "Sig","Not Sig")
data$label <- factor(data$label, levels=c("Sig","Not Sig"))
data_sig <- data %>% filter(label == "Sig")
data_selected <- data.frame()
i <- 0
for(cid in unique(data_sig$cluster.info)){
  data_up <- data_sig[data_sig$cluster.info == cid & data_sig$avg_log2FC > 0,]
  data_down <- data_sig[data_sig$cluster.info == cid & data_sig$avg_log2FC < 0,]
  
  data_top <- data_up %>% top_n(top_num,avg_log2FC)
  data_bottom <- data_down %>% top_n(top_num,-avg_log2FC)
  data_c <- rbind(data_top,data_bottom)
  i <- i + 1
  if(i == 1){
    data_selected <- data_c
  }else{
    data_selected <- rbind(data_selected,data_c)
  }
}


p <- ggplot()+
  geom_jitter(data = data,
              aes(x = cluster.info, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4,
              alpha = 0.5) +
  geom_text_repel(
    data=data_selected,
    aes(x=cluster.info,y=avg_log2FC,label=genes),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"),
    max.overlaps=50
  ) +
theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))

pdf("SC_DEGs_jitter_plots.pdf",width = 15,height=8)
print(p)
dev.off()

