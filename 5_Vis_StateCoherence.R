library(MetBrewer)
library(circlize)
library(ggplot2)
library(ggpubr)
library(colorspace)


# Uploads -----------------------------------------------------------------

sample_ls <- (read.delim("github/samples.txt", header = FALSE))$V1

gen_clusters <- as.character(unique(unlist(sapply(c(1:length(sample_ls)), function(i){
  mp_assign <- readRDS(paste("/mp_assign_full_generalized_v2/", sample_ls[i], ".rds", sep = ""))
  return(unique(mp_assign$spot_type_meta_new))
}))))

idh_type_grade <- readRDS("extend_metadata.rds")

type_grade_anno <- data.frame(sample = sample_ls,
                              type = idh_type_grade$type[match(sample_ls,idh_type_grade$Sample)],
                              tumor_grade = as.character(idh_type_grade$TumorGrade[match(sample_ls,idh_type_grade$Sample)]),
                              sample_grade = as.character(idh_type_grade$SectionGrade[match(sample_ls,idh_type_grade$Sample)]),
                              categorical_grade = idh_type_grade$CategoricalGrade[match(sample_ls,idh_type_grade$Sample)],
                              junction = idh_type_grade$junction[match(sample_ls,idh_type_grade$Sample)])


# State coherence  ------------------------------------------------------

rand_num <- 100

all_state_coh <- sapply(c(1:length(sample_ls)), function(i){
  
  print(sample_ls[i])
  
  # load data
  spots_positions <- read.csv(paste("github/data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("/mp_assign_full_generalized_v2/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  # abundance
  
  programs_comp <- sample_programs_composition(spots_clusters,gen_clusters)
  
  # neighbors tables
  
  neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
  
  rand_neighbors_table <- lapply(c(1:rand_num), function(i){
    new_pos <- sample(spots_positions$V1[spots_positions$V2 ==1], length(spots_positions$V1[spots_positions$V2 ==1]), replace = FALSE)
    pos_table <- spots_positions
    pos_table$V1[pos_table$V2 ==1] <- new_pos
    #pos_table$V2 <- spots_positions[new_pos, "V2"]
    
    neighbors_table <- neighbors_table_func(pos_table,spots_clusters)
    return(neighbors_table)
  })
  
  
  # spatial coherence
  
  programs_spatial_score <- sapply(sort(gen_clusters), function(cluster){
    if (!(cluster %in% spots_clusters$spot_type)) {
      prog_score <- NaN
    } else {
      program_neighbors_table = neighbors_table[row.names(neighbors_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster],]
      obs <- obs_program_spatial_score(program_neighbors_table, cluster)
      one_spatial_score <- one_val(dim(program_neighbors_table)[1])
      zero_spatial_score <- zero_val(rand_neighbors_table, spots_clusters, cluster)
      if (obs>one_spatial_score){obs <- one_spatial_score}
      if (obs<zero_spatial_score){obs <- zero_spatial_score}
      
      prog_score <- (obs - zero_spatial_score)/(one_spatial_score - zero_spatial_score)
    }
    return(prog_score)
  })
  
})

spatial_long <- readRDS("/results/spatial_coherence/spatial_long_genv2.rds")
# Figures -----------------------------------------------------------------

# Figure 3A

spatial_long$type1 <- ifelse(spatial_long$type == "GBM", "GBM", 
                             ifelse(spatial_long$type %in% c("Oligo2","Oligo3"), "Oligo", "Astro"))

my_comparisons <- list( c("Oligo", "Astro"), c("Oligo", "GBM"), c("Astro", "GBM"))
ggboxplot(spatial_long, x = "type1", y = "spatial_score")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1)

# Figure 3E


my_comparisons <- list( c("low", "med"), c("med", "high"), c("high", "gbm") )
ggboxplot(spatial_long, x = "cat_grade", y = "spatial_score")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5)

# Figure S3B 

samp_colors <- as.factor(ifelse(levels(spatial_long$sample) %in% type_grade_anno$sample[type_grade_anno$categorical_grade == "low"], "low",
                                ifelse(levels(spatial_long$sample) %in% type_grade_anno$sample[type_grade_anno$categorical_grade == "med"], "med",
                                       ifelse(levels(spatial_long$sample) %in% type_grade_anno$sample[type_grade_anno$categorical_grade == "high"], "high",
                                              ifelse(levels(spatial_long$sample) %in% type_grade_anno$sample[type_grade_anno$categorical_grade == "gbm"], "gbm",NA)))))


levels(samp_colors) <- c("#ee9b43","#1b7975","#78adb7","#4c9a77")
meta_pal<-c("Neuron"="#661f66","Oligo"="#C993A2","Vasc"="#762b35","Immune"="#1b7975", "Cancer"="#ee9b43","Synaptic"="#005187")

ggplot(spatial_long, aes(x=factor(sample), y=spatial_score, fill = metaprogram)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size= 12, angle = 45, vjust = 0.5, color = as.character(samp_colors)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  ylab("spatial cohernece score") +
  xlab("sample") +
  ggtitle("Spatial Coherence by samples") +
  scale_fill_manual(values = meta_pal)



# Functions ---------------------------------------------------------------

sample_programs_composition <- function(spots_clusters, gen_clusters){
  composition <- table(spots_clusters$spot_type)
  old_clusters <- names(composition)
  add_clusters <- gen_clusters[!gen_clusters %in% old_clusters]
  sapply(add_clusters,function(clust){
    composition <<- c(composition, clust = 0)
  })
  
  names(composition) <- c(old_clusters, add_clusters)
  final_composition <- composition[sort(names(composition))]/sum(composition)
  return(final_composition)
}


neighbors_table_func <- function(spots_positions,spots_clusters){
  neighbors_table <- sapply(spots_clusters$barcodes, function(spot){
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    
    if (spots_col == 0 | spots_row == 0) {
      c1 = NaN
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% spots_clusters$barcodes)){
        c1 = NaN
      } else {
        c1 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n1])
      }
    }
    
    if (spots_col == 127 | spots_row == 0) {
      c2 = NaN
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% spots_clusters$barcodes)){
        c2 = NaN
      } else {
        c2 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n2])
      }
    }
    
    if (spots_col == 0 | spots_col == 1) {
      c3 = NaN
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% spots_clusters$barcodes)){
        c3 = NaN
      } else {
        c3 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n3])
      }
    }
    
    if (spots_col == 126 | spots_col == 127) {
      c4 = NaN
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% spots_clusters$barcodes)){
        c4 = NaN
      } else {
        c4 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n4])
      }
    }
    
    if (spots_col == 0 | spots_row == 77) {
      c5 = NaN
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% spots_clusters$barcodes)){
        c5 = NaN
      } else {
        c5 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n5])
      }
    }
    
    if (spots_col == 127 | spots_row == 77) {
      c6 = NaN
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% spots_clusters$barcodes)){
        c6 = NaN
      } else {
        c6 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n6])
      }
    }
    
    
    return(c(c1,c2,c3,c4,c5,c6))
    
  })
  
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_clusters$barcodes
  
  return(neighbors_table)
}


obs_program_spatial_score <- function(program_neighbors, cluster){
  cluster_neighbors_bin <- ifelse(program_neighbors == cluster, 1, 0)
  if(is.null(dim(program_neighbors))){
    cluster_neighbors_sum <- sum(cluster_neighbors_bin)
  } else {
    cluster_neighbors_sum <- apply(cluster_neighbors_bin,1,function(rx){sum(na.omit(rx))})
  }
  obs <- mean(cluster_neighbors_sum)
  return(obs)
}


one_val <- function(spots_num){
  a <- sqrt((4*spots_num)/(6*sqrt(3)))
  oneval <- (6*spots_num-12*a-6)/spots_num
  return(oneval)
}

zero_val <- function(rand_table, spots_clusters, cluster){
  all_zeroval <- sapply(rand_table, function(neighbors_rand_table){
    program_rand_neighbors_table = neighbors_rand_table[row.names(neighbors_rand_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster],]
    rand_obs <- obs_program_spatial_score(program_rand_neighbors_table, cluster)
    return(rand_obs)
  })
  zeroval <- mean(all_zeroval)
  return(zeroval)
}

