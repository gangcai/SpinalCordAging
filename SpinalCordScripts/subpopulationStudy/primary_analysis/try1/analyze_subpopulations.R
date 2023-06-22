library(Seurat)
obj.raw <- readRDS("../../mergedAnalysis/SC_SeuratV5_CellAnnoated.rds")
#record program start time
start_time <- Sys.time()

for(cellType in levels(Idents(obj.raw))){
	obj <- subset(obj.raw,idents=cellType)

	### set the parameters for sub-population analysis
	pc.num <- 50
	resolution <- 0.3

	### normalization and scale, dimension reduction
	obj <- NormalizeData(obj)
	obj <- FindVariableFeatures(obj)
	obj <- ScaleData(obj)
	obj <- RunPCA(obj,npcs=pc.num)


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

	saveRDS(obj,paste0(cellType,"_SeuratV5.rds"))

	#### UMAP plots
	p1 <- DimPlot(
	  obj,
	  reduction = "umap",
	  group.by = c("samples","clusters"),
	  combine = T,
	  ncol=2
	)

	pdf(paste0(cellType,"_umap.samplesMerged.pdf"),width=9,height=4)
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

	pdf(paste0(cellType,"_umap_samplesSplitted.pdf"),width=8,height=6)
	print(p2)
	dev.off()


	p3 <- DimPlot(
	  obj,
	  reduction = "umap",
	  group.by = c("clusters"),
	  label=T
	)

	pdf(paste0(cellType,"_umap_showClusters.pdf"),width=5,height=4)
	print(p3)
	dev.off()

	############# Find marker genes for each cluster ###############
	obj <- JoinLayers(obj)
	markers <- FindAllMarkers(obj,
				     only.pos = TRUE,
				     assay="RNA",
				     slot="data",
				     min.pct = 0.25, logfc.threshold = 0.25)

	write.table(markers,file=paste0(cellType,"_merged_markers.tsv"),sep="\t",quote=F,row.names=F)

	############ output cell ids for each cluster ###################
	clusters=obj$seurat_clusters
	clusters=data.frame(clusters)
	write.table(clusters,file=paste0(cellType,"_cell_cluster_ids.tsv"),sep="\t",row.names=T,col.names=F,quote=F)


}
############# Print sesssion info ############
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#check elapsed time
end_time <- Sys.time()
print(capture.output(end_time-start_time))

##end of the pipeline
