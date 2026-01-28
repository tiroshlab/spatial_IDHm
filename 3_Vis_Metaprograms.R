library(reshape2)
library(ggplot2)
library(scales)
library(NMF)
library(RColorBrewer)

#####scalop functions
substri = function(x,
                   sep = "\\.|-|_",
                   pos = 1,
                   max.nchar = NULL,
                   na.rm = TRUE) {
  
  # separate strings into substrings
  strings = stringr::str_split(x, sep)
  # select relevant substrings
  strings = sapply(strings, `[`, pos, simplify = F)
  if (na.rm) {
    # remove NA positions
    strings = sapply(strings, function(s) s[!is.na(s)], simplify = F)
  }
  
  # if there is more than one substring,
  # paste selected substrings back together by the first sep delimiter
  sep = sapply(stringr::str_split(sep, "\\|"), `[`, 1)
  sep = stringr::str_replace_all(sep, "\\\\", "")
  strings = sapply(strings, paste0, collapse = sep, simplify = F)
  # convert strings from list to character vector
  strings = as.character(unlist(strings))
  # trim to keep 'max.nchar' characters per substring
  if (!is.null(max.nchar)) {
    strings = sapply(strings, function(string) {
      if (nchar(string) <= max.nchar)
        string
      else substr(string, 1, max.nchar)
    })
  }
  
  strings
}

Unlist = function(L, nested.names = FALSE) {
  if (nested.names) {
    Names = unlist(sapply(L, names), use.names = F)
  } else {
    Names = rep(names(L), lengths(L))
  }
  stats::setNames(unlist(L), Names)
}

#helper function for metaprogram generation

robust_nmf_programs <- function(nmf_programs, intra_min = 25, intra_max = 10, inter_filter=T, inter_min = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}


# Custom color palette
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

## Create list of NMF matrics where each sample is an entry
path <- "/nmf/out_comb_dir/"
sample_ls <- list.files(path)

prog_genes_ls <- readRDS(file = "/results/NMF_comb/all_nmf_wBasis_progs.RDS")

Genes_nmf_w_basis <- prog_genes_ls
nmf_programs_sig <- prog_genes_ls

# get gene programs (top 50 genes by NMF score)
nmf_programs_sig <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_filter_ccle <- robust_nmf_programs(nmf_programs_sig, intra_min = 25, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs_sig <- do.call(cbind, nmf_programs_sig)

nmf_programs_sig<-readRDS("/results/NMF_comb/programs_ls.rds")

# calculate similarity between programs
nmf_intersect <- apply(nmf_programs_sig , 2, function(x) apply(nmf_programs_sig , 2, function(y) length(intersect(x,y)))) 
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_ccle <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_ccle), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_ccle), order.dendrogram(nmf_intersect_hc_ccle)]

setwd("/results/NMF_comb/")

nmf_intersect<-readRDS("/results/NMF_comb/nmf_intersect.RDS")
nmf_programs_sig<-readRDS("/results/NMF_comb/nmf_programs_sig.RDS")

nmf_intersect_KEEP    <- nmf_intersect
nmf_programs_sig_KEEP <- nmf_programs_sig

### Parameters 

Min_intersect_initial <- 8  # the minimal intersection cutoff for defining the Founder NMF program of a cluster
Min_intersect_cluster <- 8 # the minimal intersection cuttof for adding a new NMF to the forming cluster 
Min_group_size        <- 2     # the minimal group size to consider for defining the Founder_NMF of a MP 
Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
k <- 1
Curr_cluster <- c()
MP_list      <- list()

while (Sorted_intersection[1]>Min_group_size) {   ### CHECK!
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                   <- nmf_programs_sig[,names(Sorted_intersection[1])] # initial genes are those in the first NMF. Genes_MP always has only 50 genes consisting of the current MP
  nmf_programs_sig           <- nmf_programs_sig[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs_sig))]  # remove selected NMF
  Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                <- Genes_MP  # has all genes in all NMFs in the current cluster, for newly defining Genes_MP after adding a new NMF 
  
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {   ### Define current cluster 
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- strsplit(i, "[.]")[[1]][[1]]
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history   <- c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]
    
    nmf_programs_sig      <- nmf_programs_sig[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs_sig))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  
  
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}



inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_KEEP)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters


# plot re-ordered similarity matrix heatmap     
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new,inds_new]) 

library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

Cluster_list<-readRDS("Cluster_list.rds")

MP <-  do.call(cbind, MP_list)
names(MP_list)=c("Vasc","Oligo","Neuron","cAC2","Myeloid","cUndiff","cMESHyp","cOPC1","cMES","LQ","cProlifMetab","cOXPHOS","cMESAst","cMESVasc","Synaptic","cMESInf","cNPC","cOPC2")
#saveRDS(MP_list,"combined_metaprograms_orig.rds")


####Jaccard heatmap keeping cluster order of redefined clusters
inds_sorted <- c()
for (j in 1:length(Cluster_list)){
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_KEEP)))
}
inds_sorted->inds_new
#inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters

# plot re-ordered similarity matrix heatmap     
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new,inds_new]) 


# Plot re-ordered similarity matrix heatmap
MP_plot <- reshape2::melt(nmf_intersect_KEEP[inds_new, inds_new]) 
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

## Annotate metaprograms heatmap
library(stringr)
library(dplyr)

metaprog_df <- data.frame(Programs = colnames(nmf_intersect_KEEP)[inds_new]) %>% 
  mutate(Cluster = as.factor(ifelse(Programs %in% Unlist(Cluster_list), yes = names(Unlist(Cluster_list))[match(Programs, Unlist(Cluster_list))], no = "No_Cluster")),
         Class = substri(Programs, pos = 2),
         Location = NA,
         Sample = substri(Programs, pos = 1))
levels(metaprog_df$Cluster) <- list(Vasc = "Cluster_1", Oligo = "Cluster_2", Neuron = "Cluster_3", cAC2 = "Cluster_4", Myeloid = "Cluster_5", cUndiff = "Cluster_6",
                                    cMESHyp = "Cluster_7", cOPC1 = "Cluster_8", cMES = "Cluster_9", LQ = "Cluster_10",cProlifMetab = "Cluster_11",cOXPHOS="Cluster_12", cMESAst="Cluster_13",cMESVasc="Cluster_14",Synaptic="Cluster_15",cMESInf="Cluster_16",cNPC="Cluster_17",cOPC2="Cluster_18")


##ggbar function
ggbar = function(colVar,
                 legend_title = NULL,
                 obs=NULL,
                 dir=c('h','v'),
                 cols = c('blue','red','orange','green','magenta')) {
  
  L = list(i = obs,col = colVar)
  lens = sapply(L, length, simplify = F)
  lens = lens[lens!=0]
  lens = unlist(lens)
  
  stopifnot(length(lens)>0 & length(unique(lens))==1)
  
  len = unique(unlist(lens))
  if (is.null(obs)) obs = 1:len
  if (is.null(colVar)) colVar = rep('', len)
  
  dir = match.arg(dir)
  d = data.frame(id = obs,
                 colVar=colVar,
                 stringsAsFactors=F)
  
  d = d %>% dplyr::mutate(id = factor(id, levels = unique(id)))
  
  if (dir=='h') {G = ggplot(d, aes(y=1,x=id,fill=colVar))}
  else {G = ggplot(d, aes(y=id,x=1,fill=colVar))}
  
  
  G + geom_tile() +
    scale_fill_manual(values=cols, name = legend_title, guide = guide_legend(override.aes = list(size = 10))) +
    theme_void() +
    theme(legend.position='top',
          plot.margin=margin(0,0,0,0,'cm')) 
}


# metadata
metaprog_df <- metaprog_df %>%
  mutate(IDHmut_grade = case_when(
    (Sample=="HBT277B") ~"4",(Sample=="HBT44") ~"2",(Sample=="HBT204A") ~"2",
    (Sample=="HBT314D") ~"3",(Sample=="HBT261") ~"2",(Sample=="MGH257") ~"2",
    (Sample=="HBT275") ~"3",(Sample=="HBT277A") ~"4",(Sample=="HBT312B") ~"2",
    (Sample=="HBT301B") ~"3",(Sample=="BWH24A") ~"4",(Sample=="BWH25A") ~"3",
    (Sample=="HBT301A") ~"3",(Sample=="UKF589") ~"3",(Sample=="MGH259O") ~"3",
    (Sample=="HBT314B") ~"4",(Sample=="UKF605") ~"4",(Sample=="HBT204B") ~"2",
    (Sample=="BWH35O") ~"2",(Sample=="HBT52") ~"3",(Sample=="BWH23O") ~"3",
    (Sample=="BWH28A") ~"2",(Sample=="UKF789") ~"NA",(Sample=="UKF296") ~"NA",
    (Sample=="ZH8812bulk") ~"NA",(Sample=="UKF259") ~"NA",(Sample=="UKF266") ~"NA",
    (Sample=="UKF269") ~"NA",(Sample=="UKF260") ~"NA",(Sample=="UKF334") ~"NA",
    (Sample=="UKF275") ~"NA",(Sample=="UKF248") ~"NA",(Sample=="UKF304") ~"NA",
    (Sample=="ZH8811Bbulk") ~"NA",(Sample=="ZH1019T1") ~"NA",(Sample=="ZH881T1") ~"NA",
    (Sample=="ZH1019inf") ~"NA",(Sample=="ZH8811Abulk") ~"NA",(Sample=="UKF255") ~"NA",
    (Sample=="ZH916bulk") ~"NA",(Sample=="ZH1007nec") ~"NA",(Sample=="UKF251") ~"NA",
    (Sample=="ZH916T1") ~"NA",(Sample=="UKF313") ~"NA",(Sample=="ZH881inf") ~"NA",
    (Sample=="ZH916inf") ~"NA",(Sample=="UKF243") ~"NA", TRUE ~ "NA"
  )) %>%
  mutate(grade = case_when(
    (Sample=="HBT277B") ~"4",(Sample=="HBT44") ~"2",(Sample=="HBT204A") ~"2",
    (Sample=="HBT314D") ~"3",(Sample=="HBT261") ~"2",(Sample=="MGH257") ~"2",
    (Sample=="HBT275") ~"3",(Sample=="HBT277A") ~"4",(Sample=="HBT312B") ~"2",
    (Sample=="HBT301B") ~"3",(Sample=="BWH24A") ~"4",(Sample=="BWH25A") ~"3",
    (Sample=="HBT301A") ~"3",(Sample=="UKF589") ~"3",(Sample=="MGH259O") ~"3",
    (Sample=="HBT314B") ~"4",(Sample=="UKF605") ~"4",(Sample=="HBT204B") ~"2",
    (Sample=="BWH35O") ~"2",(Sample=="HBT52") ~"3",(Sample=="BWH23O") ~"3",
    (Sample=="BWH28A") ~"2",(Sample=="UKF789") ~"NA",(Sample=="UKF296") ~"4",
    (Sample=="ZH8812bulk") ~"4",(Sample=="UKF259") ~"4",(Sample=="UKF266") ~"4",
    (Sample=="UKF269") ~"4",(Sample=="UKF260") ~"4",(Sample=="UKF334") ~"4",
    (Sample=="UKF275") ~"4",(Sample=="UKF248") ~"4",(Sample=="UKF304") ~"4",
    (Sample=="ZH8811Bbulk") ~"4",(Sample=="ZH1019T1") ~"4",(Sample=="ZH881T1") ~"4",
    (Sample=="ZH1019inf") ~"4",(Sample=="ZH8811Abulk") ~"4",(Sample=="UKF255") ~"4",
    (Sample=="ZH916bulk") ~"4",(Sample=="ZH1007nec") ~"4",(Sample=="UKF251") ~"4",
    (Sample=="ZH916T1") ~"4",(Sample=="UKF313") ~"4",(Sample=="ZH881inf") ~"4",
    (Sample=="ZH916inf") ~"4",(Sample=="UKF243") ~"4", TRUE ~ "NA"
  )) %>%
  mutate(type = case_when(
    (Sample=="HBT277B") ~"IDH-A",(Sample=="HBT44") ~"IDH-O",(Sample=="HBT204A") ~"IDH-A",
    (Sample=="HBT314D") ~"IDH-A",(Sample=="HBT261") ~"IDH-A",(Sample=="MGH257") ~"IDH-O",
    (Sample=="HBT275") ~"IDH-A",(Sample=="HBT277A") ~"IDH-A",(Sample=="HBT312B") ~"IDH-O",
    (Sample=="HBT301B") ~"IDH-A",(Sample=="BWH24A") ~"IDH-A",(Sample=="BWH25A") ~"IDH-A",
    (Sample=="HBT301A") ~"IDH-A",(Sample=="UKF589") ~"IDH-A",(Sample=="MGH259O") ~"IDH-O",
    (Sample=="HBT314B") ~"IDH-A",(Sample=="UKF605") ~"IDH-A",(Sample=="HBT204B") ~"IDH-A",
    (Sample=="BWH35O") ~"IDH-O",(Sample=="HBT52") ~"IDH-A",(Sample=="BWH23O") ~"IDH-O",
    (Sample=="BWH28A") ~"IDH-A",(Sample=="UKF789") ~"NA",
    (Sample %in% c("UKF296","ZH8812bulk","UKF259","UKF266","UKF269","UKF260","UKF334","UKF275",
                   "UKF248","UKF304","ZH8811Bbulk","ZH1019T1","ZH881T1","ZH1019inf","ZH8811Abulk",
                   "UKF255","ZH916bulk","ZH1007nec","UKF251","ZH916T1","UKF313","ZH881inf","ZH916inf","UKF243")) ~ "GBM",
    TRUE ~ "NA"
  ))

metaprog_df <- metaprog_df %>%
  mutate(
    cat_grade = case_when(
      IDHmut_grade == 2 ~ "low",
      type == "IDH-A" & IDHmut_grade == 3 ~ "mid",
      type == "IDH-O" & IDHmut_grade == 3 ~ "high",
      type == "IDH-A" & IDHmut_grade == 4 ~ "high",
      type == "GBM" ~ "GBM",
      TRUE ~ NA_character_   # fallback if none of the conditions match
    )
  )

meta_pal <- c(
  "Oligo" = "#6e948c",
  "Neuron" = "#251714",
  "cOPC1" = "#d9636c",
  "cOPC2" = "#ff9898",
  "cOXPHOS" = "#62205f",
  "cUndiff" = "#b0799a",
  "cNPC" = "#163274",
  "cProlifMetab" = "#7b89bb",
  "cAC1" = "#ef7923",
  "cAC2"="#e74b47",
  "cMESHyp" = "#FFB224",
  "cMESAst" = "#99610a",
  "cMES" = "#FFF05A",
  "cMESVasc" = "#E1BA50",
  "cMESInf" = "#FFEC9D",
  "Vasc" = "#CF1C90",
  "Myeloid" = "#5e9432",
  "Synaptic" = "#244422")

annotate_clusters <- ggbar(metaprog_df$Cluster, dir = "h", cols =meta_pal, legend_title = "Cluster") + theme(legend.direction = "horizontal")
#annotate_grade <- ggbar(metaprog_df$grade, dir = "h", cols = c("#CAB2D6","#8DD3C7","#80B1D3","darkgrey"), legend_title = "Class") + theme(legend.direction = "horizontal")
annotate_cat_grade<-ggbar(metaprog_df$cat_grade, dir = "h", cols = c("#d7aca1","#ddc000","#79ad41","#34b6c6"), legend_title = "Class") + theme(legend.direction = "horizontal")
annotate_type <- ggbar(metaprog_df$type, dir = "h", cols = c("#dd5129","#0f7ba2","#43b284","darkgrey"), legend_title = "cancer type") + theme(legend.direction = "horizontal")
annotate_samples <- ggbar(metaprog_df$Sample, dir = "h", cols = met.brewer("Signac", 50), legend_title = "Sample") + theme(legend.direction = "horizontal")

library(patchwork)
MP_plot <- annotate_cat_grade + theme(legend.position = "none") + 
  annotate_samples + theme(legend.position = "none") + 
  theme(legend.text = element_text(size = 24))+
  annotate_clusters + theme(legend.position = "none") + 
  theme(legend.text = element_text(size = 24))+
  annotate_type + theme(legend.position = "none") + 
  theme(legend.text = element_text(size = 24))+
  p1 + plot_layout(nrow = 6, heights = c(0.05, 0.05, 0.05,0.05, 1), guides = "collect")

library(grid)
lej1 <- cowplot::get_legend(annotate_clusters)
lej2 <- cowplot::get_legend(annotate_samples)
#lej3 <- cowplot::get_legend(annotate_grade)
lej4 <- cowplot::get_legend(annotate_type)
lej5 <- cowplot::get_legend(annotate_cat_grade)
grid.newpage()
grid.draw(lej1)
grid.newpage()
grid.draw(lej2)
#grid.newpage()
#grid.draw(lej3)
grid.newpage()
grid.draw(lej4)
grid.newpage()
grid.draw(lej5)

