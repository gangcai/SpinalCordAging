library(dplyr)
library(ggplot2)
##### loading the scRNA-seq data with cellType information and count matrix #########
#obj <- readRDS("../mouse_spinal_cord_c3_SeuratAnalysis.rds")  # a list with Seurat meta_data and raw counts

#raw_counts <- data_raw$counts
#meta_data <- as_tibble(obj@meta.data)
#meta_data$group <- gsub("[0-9]","",meta_data$samples)  # prepare group information based on your sample names, change needed for your samples
#meta_data$group <- gsub("-M","",meta_data$group)
#col_Data  <- meta_data %>% select(samples, group) %>% distinct()

meta_data <- read.table("../SC_MetaData.tsv",header=T,sep="\t")

cellTypes <- unique(meta_data$cellType)
cellTypes <- as.character(cellTypes)

samples <- unique(meta_data$samples)
groups <- unique(meta_data$group)


#get number of cells for each sample
sample_cellnum_list <- table(meta_data$samples)

######### get number of cells for each sample in each cellType ############
cell_num_list <- list()

for(cid in cellTypes){
  data_f <- sapply(samples, function(x){
    meta_data_f <- meta_data %>%
      filter(cellType == cid, samples == x)
    cells_f <- meta_data_f$cell
    cell_num <- length(cells_f)
    
    sample_total_cell <- sample_cellnum_list[[x]]
    
    cell_per <- cell_num/sample_total_cell
    
    
    group_id <- gsub("[0-9]","",x)
    group_id <- gsub("-M","",group_id)
    result <- c(cid, group_id, cell_num,cell_per)
    names(result) <- c("cellType_id","group_id","cell_num","cell_percentage")
    return(result)
  })
  data_f <- data.frame(t(data_f))
  cell_num_list[[cid]] <- data_f
}

######### t-test ###########
for(cid in cellTypes){
  
   df <- cell_num_list[[cid]]
   num_1 <- as.numeric(df[df$group_id == "SC-O",]$cell_percentage)
   num_2 <- as.numeric(df[df$group_id == "SC-Y",]$cell_percentage)
   fc <- log2(mean(num_1)/mean(num_2))
   p <- t.test(num_1,num_2)$p.value
   cell_num_list[[cid]]$log2_FC <- fc
   cell_num_list[[cid]]$p_value <- p
}



cell_num_tbl <- bind_rows(cell_num_list)
write.table(cell_num_tbl, file = "SC_Cell_Proportion_Changes.tsv",sep="\t",row.names = F, col.names=T, quote=F)


######## Plot the results #############
# reload data 
#cell_num_tbl <- read.table("SC_Cell_Proportion_Changes.tsv",header=T,sep="\t")
cellTypes <- unique(cell_num_tbl$cellType_id)
cellTypes <- as.character(cellTypes)
cellTypes <- cellTypes[order(cellTypes)]
cell_num_tbl$cell_percentage <- 100*as.numeric(cell_num_tbl$cell_percentage)
cell_num_tbl$cellType_id <- as.character(cell_num_tbl$cellType_id)
cell_num_tbl$cellType_id <- factor(cell_num_tbl$cellType_id, levels=as.character(cellTypes))

p <- ggplot(cell_num_tbl, aes(x = cellType_id, y = cell_percentage, fill=group_id)) +
     geom_boxplot() +
     xlab("Cell Types") +
     ylab("Cell Percentage")
pdf("SC_Cell_Proprotion_Changes_Boxplot.pdf", width=10,height=4)
print(p)
dev.off()
