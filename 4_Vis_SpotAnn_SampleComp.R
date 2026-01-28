library(scalop)
library(dplyr)
library(tidyverse)
library(tibble)
library(reshape2)
library(tidyr)
library(data.table)
library(tidytext)


####Annotate spots by dominant metaprogram###########
metaprograms_gene_list<-readRDS("/spatial_all_glioma_metaprograms.rds")


# generate normalized exp matrices
file_paths <- list.files(path = "/mats/counts_mats/", pattern = "\\.rds$", full.names = TRUE)
sample_ls  <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))
samples_names <- sample_ls

per_sample_mat <- lapply(file_paths, readRDS)

for (i in seq_along(per_sample_mat)) {
  m <- as.matrix(per_sample_mat[[i]])
  
  # remove MT/RPL/RPS genes
  m <- m[-grep("^MT-|^RPL|^RPS", rownames(m)), ]
  
  # drop zero-sum columns (spots)
  if (min(colSums(m)) == 0) {
    m <- m[, colSums(m) != 0, drop = FALSE]
  }
  
  # CPM + log transform
  scaling_factor <- 1000000 / colSums(m)
  m_CPM  <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  m_loged <- log2(1 + (m_CPM / 10))
  
  # remove genes with zero variance
  var_filter <- apply(m_loged, 1, var)
  m_proc <- m_loged[var_filter != 0, , drop = FALSE]
  
  # filter lowly expressed genes
  exp_genes <- rownames(m_proc)[rowMeans(m_proc) > 0.4]
  m_proc <- m_proc[exp_genes, , drop = FALSE]
  
  # store back
  per_sample_mat[[i]] <- m_proc
  names(per_sample_mat)[i] <- sample_ls[[i]]
  
  rm(m, m_CPM, m_loged, var_filter, exp_genes, m_proc, scaling_factor)
}

# generate score matrices (one per sample)
score_mat <- lapply(seq_along(per_sample_mat), function(i) {
  m_proc <- per_sample_mat[[i]]
  
  signatures <- scalop::sigScores(
    m_proc,
    metaprograms_gene_list,
    expr.center = TRUE,
    conserved.genes = 0.4,
    expr.nbin = 20,
    expr.binsize = 75
  )
  
  spot_scores <- data.frame(spot_names = rownames(signatures))
  spot_scores$Oligo          <- signatures$Oligo
  spot_scores$Neuron         <- signatures$Neuron
  spot_scores$cOPC1          <- signatures$cOPC1
  spot_scores$cMES           <- signatures$cMES
  spot_scores$cProlifMetab   <- signatures$cProlifMetab
  spot_scores$cOXPHOS        <- signatures$cOXPHOS
  spot_scores$cMESAst        <- signatures$cMESAst
  spot_scores$cMESVasc       <- signatures$cMESVasc
  spot_scores$Synaptic       <- signatures$Synaptic
  spot_scores$cMESInf        <- signatures$cMESInf
  spot_scores$cNPC           <- signatures$cNPC
  spot_scores$cOPC2          <- signatures$cOPC2
  spot_scores$VascAng        <- signatures$VascAng
  spot_scores$VascBBB        <- signatures$VascBBB
  spot_scores$GAM            <- signatures$GAM
  spot_scores$MacScav        <- signatures$MacScav
  spot_scores$InfMg          <- signatures$InfMg
  spot_scores$cAC1           <- signatures$cAC1
  spot_scores$cAC2           <- signatures$cAC2
  spot_scores$cMESHypLipid   <- signatures$cMESHypLipid
  spot_scores$cMESHypStress  <- signatures$cMESHypStress
  spot_scores$cUndiff1       <- signatures$cUndiff1
  spot_scores$cUndiff2       <- signatures$cUndiff2
  
  spot_scores
})

names(score_mat) <- sample_ls

# save score dfs and mp assignments per sample
for (i in seq_along(score_mat)) {
  score_df <- as.data.frame(score_mat[[i]])
  
  # save score matrix
  saveRDS(
    score_df,
    file.path("/mats/score_mats_gbm_idh_comb_generalized/", paste0(samples_names[i], ".rds"))
  )
  
  # spot assignment
  score_df <- column_to_rownames(score_df, "spot_names")
  maxcol_meta <- maxcol_strict(score_df)
  maxcol_meta <- stack(maxcol_meta)
  setnames(maxcol_meta, 2, "spot_type_meta_new")
  setnames(maxcol_meta, 1, "SpotID")
  
 
  maxcol_meta$spot_type_meta_new <- as.character(maxcol_meta$spot_type_meta_new)
  
  # metaprograms with highly overlapping gene expression collapsed to one program
  maxcol_meta$spot_type_meta_new[maxcol_meta$spot_type_meta_new %in% c("cAC1", "cAC2")] <- "cAC"
  maxcol_meta$spot_type_meta_new[maxcol_meta$spot_type_meta_new %in% c("cUndiff1", "cUndiff2")] <- "cUndiff"
  maxcol_meta$spot_type_meta_new[maxcol_meta$spot_type_meta_new %in% c("cMESHypStress", "cMESHypLipid")] <- "cMESHyp"
  maxcol_meta$spot_type_meta_new[maxcol_meta$spot_type_meta_new %in% c("InfMg", "GAM", "MacScav")] <- "Myeloid"
  maxcol_meta$spot_type_meta_new[maxcol_meta$spot_type_meta_new %in% c("VascBBB", "VascAng")] <- "Vasc"
  
  maxcol_meta$spot_type_meta_new <- factor(maxcol_meta$spot_type_meta_new)
  
  # save mp assignment
  saveRDS(
    maxcol_meta,
    file.path("/mp_assign_full_generalized_v2/", paste0(samples_names[i], ".rds"))
  )
  
  rm(score_df, maxcol_meta)
}


############sample composition by metaprograms#############

sample_ls <- (read.delim("/github/samples.txt", header = FALSE))$V1
sample_ls<-sample_ls[-c(38)] #remove UKF789 - low quality, mostly normal brain
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

gen_clusters <- as.character(unique(unlist(sapply(c(1:length(sample_ls)), function(i){
  mp_assign <- readRDS(paste("/mp_assign_full_generalized_v2/", sample_ls[i], ".rds", sep = ""))
  return(unique(mp_assign$spot_type_meta_new))
}))))

sapply(c(1:length(sample_ls)), function(i){
  
  print(sample_ls[i])
  
  # load data
  spots_clusters <- readRDS(paste("/mp_assign_full_generalized_v2/", sample_ls[i], ".rds", sep = ""))
  #spots_clusters<-rownames_to_column(spots_clusters,"barcodes")
  #select relevant column
  #spots_clusters<-spots_clusters[,c(1,20)]
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  # abundance 
  
  programs_comp <- sample_programs_composition(spots_clusters,gen_clusters)
  saveRDS(programs_comp, paste("/abund_gbm_idh_comb_generalized/", sample_ls[i],"_abund.rds", sep = ""))
})

all_comp <- list.files("/abund_gbm_idh_comb_generalized/")

compositions <- sapply(all_comp, function(comp){
  samp_comp <- readRDS(paste("/abund_gbm_idh_comb_generalized/", comp, sep = ""))
  return(samp_comp)
})


colnames(compositions) <- as.character(sapply(colnames(compositions), function(x){substr(x, 1, nchar(x)-10)}))
mp_sd <- apply(compositions,1,sd)
comp <- as.data.frame(compositions)
comp$metaprogram <- row.names(comp)
comp_long <- melt(setDT(comp), id.vars = c("metaprogram"), variable.name = "Sample", value.name = "abundance")
comp_long$metaprogram <- factor(comp_long$metaprogram, levels = unique(comp_long$metaprogram))

#metadata table v2 - sample-level annotations ####dec25 - most updated
comp_long<-comp_long %>%
  mutate(SectionGrade = case_when(
    (Sample=="HBT277B") ~"4",
    (Sample=="HBT44") ~"2",
    (Sample=="HBT204A") ~"2",
    (Sample=="HBT314D") ~"3",
    (Sample=="HBT261") ~"2",
    (Sample=="MGH257") ~"2",
    (Sample=="HBT275") ~"3",
    (Sample=="HBT277A") ~"4",
    (Sample=="HBT312B") ~"2",
    (Sample=="HBT301B") ~"3",
    (Sample=="BWH24A") ~"4",
    (Sample=="BWH25A") ~"3",
    (Sample=="HBT301A") ~"3",
    (Sample=="UKF589") ~"3",
    (Sample=="MGH259O") ~"3",
    (Sample=="HBT314B") ~"4",
    (Sample=="UKF605") ~"4",
    (Sample=="HBT204B") ~"2",
    (Sample=="BWH35O") ~"2",
    (Sample=="HBT52") ~"3",
    (Sample=="BWH23O") ~"3",
    (Sample=="BWH28A") ~"2",
    (Sample=="UKF703") ~"2",
    (Sample=="UKF243") ~"GBM",
    (Sample=="UKF248") ~"GBM",
    (Sample=="UKF251") ~"GBM",
    (Sample=="UKF255") ~"GBM",
    (Sample=="UKF259") ~"GBM",
    (Sample=="UKF260") ~"GBM",
    (Sample=="UKF266") ~"GBM",
    (Sample=="UKF269") ~"GBM",
    (Sample=="UKF275") ~"GBM",
    (Sample=="UKF296") ~"GBM",
    (Sample=="UKF304") ~"GBM",
    (Sample=="UKF313") ~"GBM",
    (Sample=="UKF334") ~"GBM",
    (Sample=="ZH1007inf") ~"GBM",
    (Sample=="ZH1007nec") ~"GBM",
    (Sample=="ZH1019inf") ~"GBM",
    (Sample=="ZH1019T1") ~"GBM",
    (Sample=="ZH8811Abulk") ~"GBM",
    (Sample=="ZH8811Bbulk") ~"GBM",
    (Sample=="ZH8812bulk") ~"GBM",
    (Sample=="ZH881inf") ~"GBM",
    (Sample=="ZH881T1") ~"GBM",
    (Sample=="ZH916bulk") ~"GBM",
    (Sample=="ZH916inf") ~"GBM",
    (Sample=="ZH916T1") ~"GBM",
    (Sample=="MGH258") ~"GBM"))

###by tumor grade
comp_long<-comp_long %>%
  mutate(type = case_when(
    (Sample=="HBT277B") ~"IDH-A",
    (Sample=="HBT44") ~"IDH-O",
    (Sample=="HBT204A") ~"IDH-A",
    (Sample=="HBT314D") ~"IDH-A",
    (Sample=="HBT261") ~"IDH-A",
    (Sample=="MGH257") ~"IDH-O",
    (Sample=="HBT275") ~"IDH-A",
    (Sample=="HBT277A") ~"IDH-A",
    (Sample=="HBT312B") ~"IDH-O",
    (Sample=="HBT301B") ~"IDH-A",
    (Sample=="BWH24A") ~"IDH-A",
    (Sample=="BWH25A") ~"IDH-A",
    (Sample=="HBT301A") ~"IDH-A",
    (Sample=="UKF589") ~"IDH-A",
    (Sample=="MGH259O") ~"IDH-O",
    (Sample=="HBT314B") ~"IDH-A",
    (Sample=="UKF605") ~"IDH-A",
    (Sample=="HBT204B") ~"IDH-A",
    (Sample=="BWH35O") ~"IDH-O",
    (Sample=="HBT52") ~"IDH-A",
    (Sample=="BWH23O") ~"IDH-O",
    (Sample=="BWH28A") ~"IDH-A",
    (Sample=="UKF703") ~"IDH-A",
    (Sample=="UKF243") ~"IDH-WT",
    (Sample=="UKF248") ~"IDH-WT",
    (Sample=="UKF251") ~"IDH-WT",
    (Sample=="UKF255") ~"IDH-WT",
    (Sample=="UKF259") ~"IDH-WT",
    (Sample=="UKF260") ~"IDH-WT",
    (Sample=="UKF266") ~"IDH-WT",
    (Sample=="UKF269") ~"IDH-WT",
    (Sample=="UKF275") ~"IDH-WT",
    (Sample=="UKF296") ~"IDH-WT",
    (Sample=="UKF304") ~"IDH-WT",
    (Sample=="UKF313") ~"IDH-WT",
    (Sample=="UKF334") ~"IDH-WT",
    (Sample=="ZH1007inf") ~"IDH-WT",
    (Sample=="ZH1007nec") ~"IDH-WT",
    (Sample=="ZH1019inf") ~"IDH-WT",
    (Sample=="ZH1019T1") ~"IDH-WT",
    (Sample=="ZH8811Abulk") ~"IDH-WT",
    (Sample=="ZH8811Bbulk") ~"IDH-WT",
    (Sample=="ZH8812bulk") ~"IDH-WT",
    (Sample=="ZH881inf") ~"IDH-WT",
    (Sample=="ZH881T1") ~"IDH-WT",
    (Sample=="ZH916bulk") ~"IDH-WT",
    (Sample=="ZH916inf") ~"IDH-WT",
    (Sample=="ZH916T1") ~"IDH-WT",
    (Sample=="MGH258") ~"IDH-WT"))

###by tumor grade
comp_long<-comp_long%>%
  mutate(TumorGrade = case_when(
    (Sample=="HBT277B") ~"4",
    (Sample=="HBT44") ~"2",
    (Sample=="HBT204A") ~"2",
    (Sample=="HBT314D") ~"3",
    (Sample=="HBT261") ~"2",
    (Sample=="MGH257") ~"3",
    (Sample=="HBT275") ~"3",
    (Sample=="HBT277A") ~"4",
    (Sample=="HBT312B") ~"3",
    (Sample=="HBT301B") ~"3",
    (Sample=="BWH24A") ~"4",
    (Sample=="BWH25A") ~"3",
    (Sample=="HBT301A") ~"3",
    (Sample=="UKF589") ~"3",
    (Sample=="MGH259O") ~"3",
    (Sample=="HBT314B") ~"3",
    (Sample=="UKF605") ~"4",
    (Sample=="HBT204B") ~"2",
    (Sample=="BWH35O") ~"2",
    (Sample=="HBT52") ~"3",
    (Sample=="BWH23O") ~"3",
    (Sample=="BWH28A") ~"2",
    (Sample=="UKF703") ~"2",
    (Sample=="UKF243") ~"GBM",
    (Sample=="UKF248") ~"GBM",
    (Sample=="UKF251") ~"GBM",
    (Sample=="UKF255") ~"GBM",
    (Sample=="UKF259") ~"GBM",
    (Sample=="UKF260") ~"GBM",
    (Sample=="UKF266") ~"GBM",
    (Sample=="UKF269") ~"GBM",
    (Sample=="UKF275") ~"GBM",
    (Sample=="UKF296") ~"GBM",
    (Sample=="UKF304") ~"GBM",
    (Sample=="UKF313") ~"GBM",
    (Sample=="UKF334") ~"GBM",
    (Sample=="ZH1007inf") ~"GBM",
    (Sample=="ZH1007nec") ~"GBM",
    (Sample=="ZH1019inf") ~"GBM",
    (Sample=="ZH1019T1") ~"GBM",
    (Sample=="ZH8811Abulk") ~"GBM",
    (Sample=="ZH8811Bbulk") ~"GBM",
    (Sample=="ZH8812bulk") ~"GBM",
    (Sample=="ZH881inf") ~"GBM",
    (Sample=="ZH881T1") ~"GBM",
    (Sample=="ZH916bulk") ~"GBM",
    (Sample=="ZH916inf") ~"GBM",
    (Sample=="ZH916T1") ~"GBM",
    (Sample=="MGH258") ~"GBM"))

#by categorical grade
comp_long<- comp_long %>%
  mutate(CategoricalGrade = case_when(
    (Sample=="HBT277B") ~"high",
    (Sample=="HBT44") ~"low",
    (Sample=="HBT204A") ~"low",
    (Sample=="HBT314D") ~"med",
    (Sample=="HBT261") ~"low",
    (Sample=="MGH257") ~"low",
    (Sample=="HBT275") ~"med",
    (Sample=="HBT277A") ~"high",
    (Sample=="HBT312B") ~"low",
    (Sample=="HBT301B") ~"med",
    (Sample=="BWH24A") ~"high",
    (Sample=="BWH25A") ~"med",
    (Sample=="HBT301A") ~"med",
    (Sample=="UKF589") ~"med",
    (Sample=="MGH259O") ~"high",
    (Sample=="HBT314B") ~"high",
    (Sample=="UKF605") ~"high",
    (Sample=="HBT204B") ~"low",
    (Sample=="BWH35O") ~"low",
    (Sample=="HBT52") ~"med",
    (Sample=="BWH23O") ~"high",
    (Sample=="BWH28A") ~"low",
    (Sample=="UKF703") ~"low",
    (Sample=="UKF243") ~"gbm",
    (Sample=="UKF248") ~"gbm",
    (Sample=="UKF251") ~"gbm",
    (Sample=="UKF255") ~"gbm",
    (Sample=="UKF259") ~"gbm",
    (Sample=="UKF260") ~"gbm",
    (Sample=="UKF266") ~"gbm",
    (Sample=="UKF269") ~"gbm",
    (Sample=="UKF275") ~"gbm",
    (Sample=="UKF296") ~"gbm",
    (Sample=="UKF304") ~"gbm",
    (Sample=="UKF313") ~"gbm",
    (Sample=="UKF334") ~"gbm",
    (Sample=="ZH1007inf") ~"gbm",
    (Sample=="ZH1007nec") ~"gbm",
    (Sample=="ZH1019inf") ~"gbm",
    (Sample=="ZH1019T1") ~"gbm",
    (Sample=="ZH8811Abulk") ~"gbm",
    (Sample=="ZH8811Bbulk") ~"gbm",
    (Sample=="ZH8812bulk") ~"gbm",
    (Sample=="ZH881inf") ~"gbm",
    (Sample=="ZH881T1") ~"gbm",
    (Sample=="ZH916bulk") ~"gbm",
    (Sample=="ZH916inf") ~"gbm",
    (Sample=="ZH916T1") ~"gbm",
    (Sample=="MGH258") ~"gbm"))

##junction type
comp_long<- comp_long %>%
  mutate(junction = case_when(
    (Sample=="HBT277B") ~"none",
    (Sample=="HBT44") ~"GGMB",
    (Sample=="HBT204A") ~"WGMJ",
    (Sample=="HBT314D") ~"WWMB",
    (Sample=="HBT261") ~"WGMJ",
    (Sample=="MGH257") ~"GGMB",
    (Sample=="HBT275") ~"none",
    (Sample=="HBT277A") ~"none",
    (Sample=="HBT312B") ~"none",
    (Sample=="HBT301B") ~"WWMB",
    (Sample=="BWH24A") ~"none",
    (Sample=="BWH25A") ~"WGMJ",
    (Sample=="HBT301A") ~"WWMB",
    (Sample=="UKF589") ~"none",
    (Sample=="MGH259O") ~"GWMJ",
    (Sample=="HBT314B") ~"none",
    (Sample=="UKF605") ~"none",
    (Sample=="HBT204B") ~"WGMJ",
    (Sample=="BWH35O") ~"GWMJ",
    (Sample=="HBT52") ~"none",
    (Sample=="BWH23O") ~"none",
    (Sample=="BWH28A") ~"none",
    (Sample=="UKF703") ~"WGMJ",
    (Sample=="UKF243") ~"gbm",
    (Sample=="UKF248") ~"gbm",
    (Sample=="UKF251") ~"gbm",
    (Sample=="UKF255") ~"gbm",
    (Sample=="UKF259") ~"gbm",
    (Sample=="UKF260") ~"gbm",
    (Sample=="UKF266") ~"gbm",
    (Sample=="UKF269") ~"gbm",
    (Sample=="UKF275") ~"gbm",
    (Sample=="UKF296") ~"gbm",
    (Sample=="UKF304") ~"gbm",
    (Sample=="UKF313") ~"gbm",
    (Sample=="UKF334") ~"gbm",
    (Sample=="ZH1007inf") ~"gbm",
    (Sample=="ZH1007nec") ~"gbm",
    (Sample=="ZH1019inf") ~"gbm",
    (Sample=="ZH1019T1") ~"gbm",
    (Sample=="ZH8811Abulk") ~"gbm",
    (Sample=="ZH8811Bbulk") ~"gbm",
    (Sample=="ZH8812bulk") ~"gbm",
    (Sample=="ZH881inf") ~"gbm",
    (Sample=="ZH881T1") ~"gbm",
    (Sample=="ZH916bulk") ~"gbm",
    (Sample=="ZH916inf") ~"gbm",
    (Sample=="ZH916T1") ~"gbm",
    (Sample=="MGH258") ~"gbm"))

saveRDS(comp_long,"/idh_gbm_comb_generalized_comp_table.rds")
