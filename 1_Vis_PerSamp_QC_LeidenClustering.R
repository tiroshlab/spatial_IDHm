setwd("/spatial_IDH/github/")
Sys.setenv(RETICULATE_PYTHON = "/.conda/envs/leiden/bin/python")

library(reticulate)
library(Seurat)
library(scalop)
library(leiden)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)
library(reshape2)

samples_names <- (read.delim("/spatial_IDH/github/samples_idh_only.txt", header = FALSE, sep = "\t"))$V1
#parameters 
complexity_filter <- 1000
mitochondrial_filter <- 20
genes_filter <- 7000
n_dim <- 20
dim_filter <- 5^(-15)
res_param <- 1
sig_th <- .005
mp_num <- 13
dim2use_list <- c(7,12,12,5,5,3,7,5,16,6,7,11,7,3,7,3,7,4,5,7,6,9,4,5) ###can set dims by running jackstraw per sample after PCA

#####Per sample QC, dimensionality reduction, Leiden clustering - generate per sample Leiden cluster gene programs
###IDH-M samples only - for GBM samples, the same approach is used - see tiroshlab/Spatial_Glioma github

leiden_clustering <- sapply(c(1:length(samples_names)), function(i){
  print(samples_names[i])
  # load spatial data
  unfilt_obj<-Load10X_Spatial(data.dir = paste("IDH_mut_data/",samples_names[i],"/outs", sep = ""),
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "detected_tissue_image.jpg",
                              filter.matrix = TRUE,
                              to.upper = FALSE
  )
  
  # filtering
  exp_obj <- subset(unfilt_obj, subset = nCount_Spatial > complexity_filter) # filter out spots with less than #complexity_filter UMIs
  exp_obj[["percent.mt"]] <- PercentageFeatureSet(exp_obj, pattern = "^MT-") # filter spots with high percent (mitochondrial_filter) mitochondrial genes
  exp_obj <- subset(exp_obj, subset = percent.mt<mitochondrial_filter)
  
  count_mat <- as.matrix(GetAssayData(exp_obj, slot = "counts"))
  #saveRDS(count_mat, paste("counts_mats_IDH_v2/", samples_names[i], "counts.rds", sep =""))
  
  # normalization and centering
  exp_obj <- NormalizeData(exp_obj, assay = "Spatial", normalization.method="LogNormalize", scale.factor=10000)
  exp_obj <- FindVariableFeatures(exp_obj, selection.method = "vst", nfeatures = genes_filter) # use only #genes_filter most variable genes
  exp_obj <- ScaleData(exp_obj)
  
  # PCA
  exp_obj <- RunPCA(exp_obj, features = VariableFeatures(object = exp_obj))
  #exp_obj <- JackStraw(exp_obj, num.replicate = 100)
  #exp_obj <- ScoreJackStraw(exp_obj, dims = 1:n_dim)
  #js_scores <- exp_obj@reductions$pca@jackstraw$overall.p.values
  #dim2use <- max(js_scores[,"PC"][js_scores[,"Score"] < dim_filter])
  # 
  dim2use <- dim2use_list[i]
  
  # clustering
  exp_obj <- FindNeighbors(exp_obj, dims = 1:dim2use)
  exp_obj <- FindClusters(exp_obj, algorithm=4, resolution = res_param) # leiden based clustering
  #print(SpatialPlot(exp_obj,label = FALSE, stroke=0, image.alpha = 0, crop = TRUE, pt.size.factor = 1.8, cols=distinct16_pal) + scale_fill_manual(values=distinct16_pal))
  
  spatial_sample_programs <- FindAllMarkers(exp_obj, only.pos = TRUE, return.thresh = sig_th)
  
  genes_num <- table(spatial_sample_programs$cluster)
  spots_clusters <- exp_obj@meta.data[,"seurat_clusters"]
  names(spots_clusters) <- row.names(exp_obj@meta.data)
  cluster_num <- table(exp_obj@meta.data$seurat_clusters)
  
  return(list(spatial_sample_programs,genes_num,spots_clusters,cluster_num)) # output: #genes per cluster, #spots per cluster, clusters genes sig, spots assignment to cluster
  
})


# cluster quality
genes_num <- c()
spots_num <- c()
sapply(c(1:length(samples_names)), function(i){
  genes_num <<- c(genes_num,as.numeric(leiden_clustering[2,i][[1]]))
  spots_num <<- c(spots_num,as.numeric(leiden_clustering[4,i][[1]]))
})
hist(genes_num, breaks = 30)
hist(spots_num, breaks = 30)

# create gene program table from Leiden clusters
genes_sig <- leiden_clustering[1,1][[1]]
genes_sig$sample <- rep(samples_names[1],dim(genes_sig)[1])

sapply(c(2:length(samples_names)),function(i){
  genes_sig_temp <- leiden_clustering[1,i][[1]]
  genes_sig_temp$sample <- rep(samples_names[i],dim(genes_sig_temp)[1])
  genes_sig <<- rbind(genes_sig,genes_sig_temp)
})

# create cluster assignment table 
spots_assign <- as.data.frame(leiden_clustering[3,1][[1]])
spots_assign$sample <- rep(samples_names[1],dim(spots_assign)[1])
row.names(spots_assign) <- paste(row.names(spots_assign), "_", samples_names[1], sep = "")
colnames(spots_assign) <- c("cluster","sample")

sapply(c(2:length(samples_names)),function(i){
  spots_assign_temp <- as.data.frame(leiden_clustering[3,i][[1]])
  spots_assign_temp$sample <- rep(samples_names[i],dim(spots_assign_temp)[1])
  row.names(spots_assign_temp) <- paste(row.names(spots_assign_temp), "_", samples_names[i], sep = "")
  colnames(spots_assign_temp) <- c("cluster","sample")
  spots_assign <<- rbind(spots_assign,spots_assign_temp)
})

# format Leiden clusters gene program list
genes_list <- lapply(c(1:length(samples_names)), function(i){
  sample_table <- leiden_clustering[1,i][[1]]
  sample_genes <- lapply(c(1:length(unique(sample_table$cluster))),function(j){
    cluster_table <- sample_table[sample_table$cluster == j,]
    genes <- cluster_table$gene[order(cluster_table$avg_log2FC, decreasing = TRUE)][1:min(50,nrow(cluster_table))]
    return(genes)
  })
  return(sample_genes)
})
genes_list <- unlist(genes_list, recursive = FALSE)
names(genes_list) <- unlist(sapply(c(1:length(samples_names)),function(i){
  paste(samples_names[i],c(1:length(leiden_clustering[4,i][[1]])), sep = "_")
}))

samples_list <- unlist(sapply(c(1:length(samples_names)),function(i){
  rep(samples_names[i],length(leiden_clustering[4,i][[1]]))
}))
saveRDS(genes_list,"all_leiden_programs.rds")












