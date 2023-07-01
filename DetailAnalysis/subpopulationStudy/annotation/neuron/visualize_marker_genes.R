library(Seurat)
obj <- readRDS("../../analysis/Neuron_SeuratV5.rds")

#gene	cellType
#Slc17a6	Excitatory
#Slc5a7	Cholinergic
#Chat	Cholinergic
#Slc6a5	Inhibitory
#Gad1	Inhibitory
#Gad2	Inhibitory

genes <- c("Slc17a6","Slc5a7","Chat","Slc6a5","Gad1","Gad2")

p <- FeaturePlot(obj, genes, reduction="umap")
pdf("Neuron_subpopu_markerGenes_FeaturePlot.pdf",width=12,height=12)
print(p)
dev.off()


p <- VlnPlot(obj, genes, ncol=1)
pdf("Neuron_subpopu_markerGenes_VlnPlot.pdf",width=6,height=12)
print(p)
dev.off()

p <- VlnPlot(obj, genes, ncol=1, pt.size=0)
pdf("Neuron_subpopu_markerGenes_VlnPlot_NoDots.pdf",width=6,height=12)
print(p)
dev.off()
