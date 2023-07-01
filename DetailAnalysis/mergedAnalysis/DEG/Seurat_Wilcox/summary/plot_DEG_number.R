library(ggplot2)
library(dplyr)

df <- read.table("DEG_Summary.tsv",header=T,sep="\t")

df_total <- df %>% group_by(cellType) %>% summarise(DEGNum = sum(DEGNum))

df_total <- df_total[order(df_total$DEGNum,decreasing=T),]

df$cellType <- factor(df$cellType, levels=df_total$cellType)
p <- ggplot(df, aes(x=cellType, y=DEGNum, fill=DEGType)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8)) +
  scale_fill_manual(values=c("darkseagreen","coral"))

pdf("SC_DEG_Summary.pdf",width=6,height=4)
print(p)
dev.off()
