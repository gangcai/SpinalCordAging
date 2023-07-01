library(Seurat)
obj <- readRDS("../primaryAnalysis/SC_SeuratV5.rds")
#cluster cellType
anno_df <- read.table("../primaryAnalysis/annotation/primary_annotation.tsv",header=T,sep="\t")
meta_data <- obj@meta.data

cellType <- sapply(meta_data$seurat_clusters,function(x){
        x <- as.character(x)
	anno_df[anno_df$cluster == x,]$cellType
})
meta_data$cellType <- unlist(cellType)

#merged and annotated idents
new_idents <- meta_data$cellType
names(new_idents) <- rownames(meta_data)
new_idents <- factor(new_idents)
Idents(obj) <- new_idents
obj@meta.data <- meta_data
saveRDS(obj,"SC_SeuratV5_CellAnnoated.rds")
write.table(meta_data,file="SC_MetaData.tsv",sep="\t",row.names=F,col.names=T,quote=F)

