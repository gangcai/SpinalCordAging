library(dplyr)
neuron_names_df <- read.table("neuron/SC_Neuron_Subpopulation_Annotation.tsv",sep="\t",header=T)
files.list <- list.files("../analysis/","cell_cluster_ids.tsv")

brief_names_list <- list("Astrocyte"="Ast",
                         "EC"="EC",
                         "Ependymal"="Epe",
                         "Meningeal"="Men",
                         "Microglia"="Mic",
                         "OL"="OL",
                         "OPC"="OPC",
                         "Pericyte"="Per",
                         "Schwann"="Sch")

i <- 0
results <- data.frame()
for(f in files.list){
  cellType <- sub("_cell_cluster_ids.tsv","",f)
  filename <- paste0("../analysis/",f)
  cluster_df <- read.table(filename,sep="\t",header=F)
  clusters <- unique(cluster_df$V2)
  if(cellType == "Neuron"){
    cid_name <- data.frame("cluster_id"=neuron_names_df$cluster_id,
                           "sub_cellType"=neuron_names_df$cell_type_name)
  }else{
    cell_names <- sapply(clusters, function(x){
      paste0(brief_names_list[[cellType]],as.numeric(as.character(x))+1)
    })
    cid_name <- data.frame("cluster_id"=clusters,
                           "sub_cellType"=cell_names)
  }
  cid_name$cellType <- cellType
  i <- i + 1
  if(i == 1){
    results <- cid_name
  }else{
    results <- rbind(results, cid_name)
  }
}

write.table(results,file = "SC_subpopulation_names.tsv",sep="\t",row.names=F,col.names=T,quote=F)
