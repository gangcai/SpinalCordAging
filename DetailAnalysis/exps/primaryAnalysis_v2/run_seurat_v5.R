### Seurat based analysis of mouse spinal cord aging (Xie lab)###

#record program start time
start_time <- Sys.time()

#loading seurat library
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(future)
plan("multicore", workers = 20)

options(future.globals.maxSize = 20000 * 1024^2) # 20G memory, Notice: change the memory usage accordingly
options(Seurat.object.assay.version = "v5") # use Seurat V5 format


### parameter settings ####
project_id <- "SC"
min.features <- 500
min.cells <- 3
pc.num <- 60         # number of PCs
resolution <- 0.4    # higher values will generate more clusters, and vice versa, recommended values for experiments: 0.05, 0.1, 0.15, 0.2, 0.25, 0.3
nCount_RNA.max <- 15000    # Notice: change this value according to your data QC results
nCount_RNA.min <- 1000
max.features <- 5000       # Notice: change this value according to your data QC results
per.mt <- 5  # maximal mitochondrial percentage

########## loading the 10X matrix files ##################
samples=c("SC-O-21M","SC-O-23M","SC-Y-4M","SC-Y-6M")    # Notice: Change the sample names to your own sample names
raw.data.merged=""

i=0
dir_base <- "/home/db/private/XieLab/aging/mouse_spinal_cord/run_celescope1.11.0/"
for(sample in samples){
        dir <- paste0(dir_base,"/",sample,"/",sample,"/05.count/",sample,"_filtered_feature_bc_matrix/")
        obj.counts <- Read10X(data.dir = dir)
        colnames=colnames(obj.counts)
        colnames(obj.counts)=paste0(sample,"_",colnames)
        i=i+1
        if(i == 1){
                raw.data.merged=obj.counts
        }else{
                raw.data.merged=cbind(raw.data.merged,obj.counts)
        }
}


all.cells=colnames(raw.data.merged)
all.samples=sapply(all.cells,function(x){strsplit(as.character(x),"_")[[1]][1]})
all.samples=as.character(all.samples)
groups <- sapply(all.samples,function(x){
  if(x == "SC-O-21M"){
    return("Old")
  }
  if(x == "SC-O-23M"){
    return("Old")
  }

  if(x=="SC-Y-4M"){
    return("Young")
  }

  if(x=="SC-Y-6M"){
    return("Young")
  }
})


metadata = data.frame("groups"=groups,"samples"=all.samples,"cells"=all.cells,row.names=all.cells)

################## Seurat basic analysis pipeline ##################
obj <- CreateSeuratObject(raw.data.merged, meta.data = metadata, project = project_id, min.cells = min.cells, min.features = min.features)

obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")


######## QC metrics #########
#-> summarize mitochondrial percentage
mt.pert <- obj@meta.data$percent.mt
mt.qt <- quantile(mt.pert,probs=seq(0,1,0.01))
mt.qt.df <- data.frame("proportion"=mt.qt, "percentile"=names(mt.qt))
mt.qt.df$percentile <- factor(mt.qt.df$percentile, levels=mt.qt.df$percentile)
p <- ggplot(mt.qt.df, aes(x=percentile, y=proportion)) + geom_bar(stat="identity") +
   theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8)) + ylab("percent.mt")+
  geom_hline(yintercept = 5,color="red")
pdf(paste0(project_id,"_percent.mt_quantile.pdf"),width=15,height=5)
print(p)
dev.off()

mt.total <- length(mt.pert)
mt.stats <- data.frame("cellPercentage"=100*c(length(mt.pert[mt.pert<2.5])/mt.total,
                         length(mt.pert[mt.pert<5])/mt.total,
                         length(mt.pert[mt.pert<7.5])/mt.total,
                         length(mt.pert[mt.pert<10])/mt.total,
                         length(mt.pert[mt.pert<20])/mt.total),
                       "percent.mt_cutoff"=factor(c("<2.5","<5","<7.5","<10","<20"),
                                                  levels=c("<2.5","<5","<7.5","<10","<20"))
                       )

p <- ggplot(mt.stats, aes(x=percent.mt_cutoff, y=cellPercentage)) + geom_bar(stat="identity",fill="lightblue") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8)) + ylab("cellPercentage") + NoLegend()
pdf(paste0(project_id,"_percent.mt_stats.pdf"),width=7,height=5)
print(p)
dev.off()


p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf(paste0(project_id, "_sctransform_QC_before.pdf"))
print(p0)
dev.off()

####### Data filtering and clustering ##########
#data filtering
obj <- subset(obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < per.mt & nCount_RNA < nCount_RNA.max &  nCount_RNA > nCount_RNA.min )


p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf(paste0(project_id, "_sctransform_QC_after.pdf"))
print(p0)
dev.off()


### split by samples
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)

### normalization and scale, dimension reduction
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj,npcs=pc.num)
obj <- FindNeighbors(obj, dims = 1:pc.num, reduction = "pca")
obj <- FindClusters(obj, resolution = resolution, cluster.name = "unintegrated_clusters")

#check un-integrated results
obj <- RunUMAP(obj, dims = 1:pc.num, reduction = "pca", reduction.name = "umap.unintegrated")
p <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("samples", "groups"))
pdf(paste0(project_id, "_umap.unintegrated.pdf"),width = 12,height=6)
print(p)
dev.off()

p <- DimPlot(obj, reduction = "umap.unintegrated", split.by = "samples")
pdf(paste0(project_id, "_umap.unintegrated_splitted.pdf"),width = 12,height=4)
print(p)
dev.off()



### intergration , slow step
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

### clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:pc.num)
obj <- FindClusters(obj, resolution = resolution, cluster.name = "clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:pc.num, reduction.name = "umap")

p1 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = c("samples","clusters"),
  combine = T,
  ncol=2
)

pdf(paste0(project_id,"_umap.samplesMerged.pdf"),width=9,height=4)
print(p1)
dev.off()


p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by=c("clusters"),
  split.by="samples",
  combine = FALSE,
  ncol=2
)

pdf(paste0(project_id,"_umap_harmony_samplesSplitted.pdf"),width=8,height=6)
print(p2)
dev.off()


p3 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = c("clusters"),
  label=T
)

pdf(paste0(project_id,"_umap_harmony_showClusters.pdf"),width=5,height=4)
print(p3)
dev.off()

############## draw for each clusters  #######################
all_idents <- Idents(obj)

p_list <- list()

clusters <- unique(all_idents)

clusters <- clusters[order(clusters)]

if(!dir.exists("show_EachCellTypeCells/")){
  dir.create("show_EachCellTypeCells/")
}

for(celltype in clusters){
  each_cells <- names(all_idents[all_idents == celltype])
  
  p <- DimPlot(obj,
               reduction = "umap",
               split.by = "samples",
               cells.highlight = each_cells) + ggtitle(paste0("c",celltype))
  pdf(paste0("show_EachCellTypeCells/",celltype,"_highlight_splitbysamples.pdf"),width=20,height=5)
  print(p)
  dev.off()
  
  p <- DimPlot(obj,
               reduction = "umap",
               cells.highlight = each_cells) + ggtitle(paste0("c",celltype))
  pdf(paste0("show_EachCellTypeCells/",celltype,"_highlight.pdf"),width=8,height=8)
  print(p)
  dev.off()
  p_list[[celltype]] <- p
}



p <- wrap_plots(p_list,ncol=4)
pdf("show_EachCellTypeCells/celltype_highlighted.pdf",width=20,height=25)
print(p)
dev.off()


#-> check umap.unintegrated
if(!dir.exists("show_EachCellTypeCells_unintegrated/")){
  dir.create("show_EachCellTypeCells_unintegrated/")
}

for(celltype in clusters){
  each_cells <- names(all_idents[all_idents == celltype])
  
  p <- DimPlot(obj,
               reduction = "umap.unintegrated",
               split.by = "samples",
               cells.highlight = each_cells) + ggtitle(paste0("c",celltype))
  pdf(paste0("show_EachCellTypeCells_unintegrated/",celltype,"_highlight_splitbysamples.pdf"),width=20,height=5)
  print(p)
  dev.off()
  
  p <- DimPlot(obj,
               reduction = "umap.unintegrated",
               cells.highlight = each_cells) + ggtitle(paste0("c",celltype))
  pdf(paste0("show_EachCellTypeCells_unintegrated/",celltype,"_highlight.pdf"),width=8,height=8)
  print(p)
  dev.off()
  p_list[[celltype]] <- p
}



p <- wrap_plots(p_list,ncol=4)
pdf("show_EachCellTypeCells_unintegrated/celltype_highlighted.pdf",width=20,height=25)
print(p)
dev.off()


######### QC metrics for each cluster after clustering ###########
p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size=0)
pdf(paste0(project_id, "_sctransform_QC_after_identssplitted.pdf"),width=12,height=8)
print(p0)
dev.off()

meta.data <- obj@meta.data

for(sample in unique(meta.data$samples)){
  sample.cells <- meta.data[meta.data$samples == sample,]$cells
  sample.obj <- subset(obj, cells = sample.cells)
  p0 <- VlnPlot(sample.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size=0)
  pdf(paste0(project_id,"_",sample, "_sctransform_QC_after_identssplitted.pdf"),width=12,height=8)
  print(p0)
  dev.off()
}


## save object
saveRDS(obj,paste0(project_id,"_SeuratV5.rds"))


############# Find marker genes for each cluster ###############

obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj,
                             only.pos = TRUE, 
                             assay="RNA",
                             slot="data",
                             min.pct = 0.25, logfc.threshold = 0.25)

top_marker_df <- markers %>% group_by(cluster) %>% top_n(5,avg_log2FC)
top_genes <- unique(top_marker_df$gene)
p <- DotPlot(obj,top_genes,assay="RNA") + theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))
pdf(paste0(project_id,"_Top_Marker_DotPlot.pdf"),height=7,width=30)
print(p)
dev.off()

write.table(markers,file=paste0(project_id,"_merged_markers.tsv"),sep="\t",quote=F,row.names=F)

############## output cell ids for each cluster ###################
clusters=obj$seurat_clusters
clusters=data.frame(clusters)
write.table(clusters,file=paste0(project_id,"_cell_cluster_ids.tsv"),sep="\t",row.names=T,col.names=F,quote=F)




################## Print sesssion info ###################
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#check elapsed time
end_time <- Sys.time()
print(capture.output(end_time-start_time))

##end of the pipeline
