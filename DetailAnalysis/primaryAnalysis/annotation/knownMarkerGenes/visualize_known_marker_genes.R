library("Seurat")
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
obj <- readRDS("../../SC_SeuratV5.rds")
DefaultAssay(obj) <- "RNA"
  
avg_exp <- AverageExpression(obj)$RNA


### loading knownMarker Genes
marker_genes_df <- read.table("knownMarkerGenes.tsv",header=T,sep="\t")
all_genes <- rownames(avg_exp)
marker_genes_df_s <- marker_genes_df[marker_genes_df$MarkerGene %in% all_genes,]

avg_exp_s <- avg_exp[marker_genes_df_s$MarkerGene,]

row_anno <- data.frame(cellType = marker_genes_df_s$CellType)
rownames(row_anno) <- marker_genes_df_s$MarkerGene

row_anno$cellType <- factor(row_anno$cellType, levels = unique(row_anno$cellType))

pdf("knownMarkerGenes_heatmap_rowscaled.pdf",height=8,width=7)
p <- pheatmap(avg_exp_s,scale="row",
              cluster_rows=F,
              annotation_row=row_anno,
              color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                        "RdBu")))(100))
print(p)
dev.off()

pdf("knownMarkerGenes_heatmap_colscaled.pdf",height=8,width=7)
p <- pheatmap(avg_exp_s,scale="column",
              cluster_rows=F,
              annotation_row=row_anno,
              color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                        "RdBu")))(100))
print(p)
dev.off()

# known marker gene dotplot
p <- DotPlot(obj, assay = "RNA", features = marker_genes_df_s$MarkerGene) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))
pdf("knownMakerGenes_Dotplot.pdf",width = 12,height=5)
print(p)
dev.off()
