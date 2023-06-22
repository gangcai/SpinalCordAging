library(Seurat)
obj <- readRDS("../../analysis/Microglia_SeuratV5.rds")
FeaturePlot(obj, c("Ctss","Ptprc","Cd3e","Msr1"), reduction="umap")
VlnPlot(obj, c("Ctss","Ptprc","Cd3e","Msr1","Cd247"))