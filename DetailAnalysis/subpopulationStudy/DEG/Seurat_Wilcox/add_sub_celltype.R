sub_celltype_df <- read.table("../../annotation/SC_subpopulation_names.tsv",header=T,
                              sep="\t")


deg_df <- read.table("SC_subpopulations_DEGs_WilCox.tsv",header=T,sep="\t")

sub_cells <- c()
for(n in c(1:nrow(deg_df))){
  celltype <- deg_df[n,"cellType"]
  cid <- deg_df[n,"cluster.info"]
  sub_celltype <- sub_celltype_df[sub_celltype_df$cluster_id == cid & sub_celltype_df$cellType == celltype,]$sub_cellType
  sub_cells <- c(sub_cells ,sub_celltype)
}

deg_df$sub_cellType <- sub_cells
write.table(deg_df,"SC_subpopulations_DEGs_WilCox_annotated.tsv",sep="\t",row.names = F, col.names=T,quote=F)