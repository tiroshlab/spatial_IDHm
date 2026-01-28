####NMF per sample#####
library(scales)
library(reshape2)
library(scales)
library(NMF)
library(MetBrewer)
library(colorspace)
library(tibble)
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(Matrix)
library(bigmemory)
library(doMC)
library(patchwork)
library(readr)

######generating input matrices for NMF
# vector of file paths
samples_names <- (read.delim("/spatial_IDH/github/samples.txt", header = FALSE, sep = "\t"))$V1
file_paths <- list.files(path ="/counts_mats/", pattern = "\\.rds", full.names = TRUE)#directory with counts mats
sample_ls <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))
per_sample_mat <- lapply(file_paths, readRDS)

for (i in seq_along(per_sample_mat)){ 
  m <- as.matrix(per_sample_mat[[i]])
  if(min(colSums(m)) == 0){m <- m[, colSums(m) != 0]}
  scaling_factor <- 1000000/colSums(m)
  m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  m_loged <- log2(1 + (m_CPM/10))
  
  # removing genes with zero variance across all cells
  var_filter <- apply(m_loged, 1, var)
  m_proc <- m_loged[var_filter != 0, ]
  # filtering out lowly expressed genes
  exp_genes <- rownames(m_proc)[(rowMeans(m_proc) > 0.4)]
  m_proc <- m_proc[exp_genes, ]
  
  # centering data gene-wise
  count_mat_cent<- m_proc - rowMeans(m_proc)
  #nmf preprocessing
  count_mat_nmf<- count_mat_cent
  count_mat_nmf[count_mat_nmf<0] <- 0 # negative values should be set to 0 in initial matrix
  
  # output to a list of gene expression profiles (GEP)
  per_sample_mat[[i]] <- count_mat_nmf
  names(per_sample_mat)[i] <- sample_ls[[i]]
  rm(m,m_loged, var_filter, exp_genes, m_proc,count_mat_cent)
  saveRDS(count_mat_nmf, paste("/nmf_dir/mat_dir/", samples_names[i], ".rds", sep =""))
}

## Extract sample names for sample list file:
samples <- list.files("/nmf_dir/mat_dir/") %>% substri(., pos = 1) %>% str_c(., collapse = "\n")
write_lines(samples, ("/nmf_dir/samples.txt"), sep = "\n")


###### NMF (per sample stand-alone local option)---------------------------------------------------------------------

#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE)

sname <- "sample"# enter sample name here 
rank_lb <- 2
rank_ub <- 10

m <- readRDS(paste0("MP/NMF/NMF_mats/", sname, ".rds"))
m <- as.matrix(m)
res <- NMF::nmf(x = m, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = list(debug=F, parallel=F, shared.memory=F, verbose=T))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)



###NMF per sample - running per sample on cluster


#inputs for running NMF (per sample): 

#txt file with sample names
#directory containing per sample gene expr matrix (centered log normalized expression matrix with all negative values changed to zeroes--non-negative), each per sample expression matrix saved as rds
sname <- as.character(args[1]) #sample name
mat_dir <- as.character(args[2]) #directory containing per sample nonnegative exp mats 
out_dir <- as.character(args[3]) #output directory
rank_lb <- as.numeric(args[4]) #lowest rank (k)
rank_ub <- as.numeric(args[5]) #highest rank (k)


print(paste0("Reading NMF matrix for sample ", sname, " from ", mat_dir))
path <- "/nmf_dir" #
m <- readRDS(paste0("/nmf_dir/mat_dir/", sname, ".rds")) 
#gene expr matrix for NMF input contains centered log normalized expression matrix with all negative value changed to zeroes (non-negative)

print(paste0("Done, matrix has ", nrow(m), " rows and ", ncol(m), " columns"))

m <- as.matrix(m)

print(paste0("Running NMF for sample ", sname, ", with the ranks- ", rank_lb, ":", rank_ub, ", started at ", Sys.time()))
#res <- NMF::nmf(x = m, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = "p")
res <- NMF::nmf(x = m, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = list(debug=F, parallel=F, shared.memory=F, verbose=T))
print(paste0("Finished at ", Sys.time()))

filename <- paste0(path, "/", out_dir, "/", sname, "_nmf_res.RDS")

print(paste0("Writing result to ", filename))

saveRDS(res, filename)

print("Done!")

###example bash script for running above R script
#!/usr/bin/env bash

# file=$1
# mat_dir=$2
# out_dir=$3
# 
# echo samples file: $file
# echo matrices directory: $mat_dir
# echo output directory: $out_dir
# 
# module load R/4.1.1.rstudio-foss-2021a
# module load blas/3.8.0
# module load gcc
# 
# echo starting...
# 
# while IFS= read -r line || [[ -n "$line" ]]; do
# sname=$line
# echo "Sending request to process sample: $sname"
# 
# bsub -q new-medium -J $sname -o $sname.o -e $sname.e -R "rusage[mem=4GB] span[hosts=1]" Rscript nmf_Rscript.R $sname $mat_dir $out_dir 2 10
# done < "$file"
# 
# echo done!
  


