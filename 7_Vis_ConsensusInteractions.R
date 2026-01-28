#######Consensus interactions########

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
})


# Helpers

# sort each cell type pair so it is consistent between analyses
sort_pair <- function(pairs) {
  sapply(pairs, function(x) paste(sort(unlist(strsplit(x, " "))), collapse = " "))
}

# rescale to [-1,1] within group 
rescale_signed_by_group <- function(df, value_col, group_col = "analysis", out_col = "pair_mean_rescaled") {
  df[[out_col]] <- NA_real_
  groups <- unique(df[[group_col]])
  for (grp in groups) {
    vals <- df[[value_col]][df[[group_col]] == grp]
    val_min <- min(vals, na.rm = TRUE)
    val_max <- max(vals, na.rm = TRUE)
    if (!is.finite(val_min) || !is.finite(val_max) || val_max == val_min) next
    
    rescaled_vals <- ifelse(
      vals < 0,
      ((vals - val_min) / (val_max - val_min)) - 1,
      ((vals - val_min) / (val_max - val_min))
    )
    df[[out_col]][df[[group_col]] == grp] <- rescaled_vals
  }
  df
}

# compute consensus lists from  "summary_df" that has:
#   rownames (pair id), type in {colocal, adjacency, proximity}, prop, mean_scaled
compute_consensus <- function(summary_df,
                              pair_col = "rownames",
                              type_col = "type",
                              prop_col = "prop",
                              score_col = "mean_scaled",
                              score_thresh = 0.35,
                              prop_thresh  = 0.20,
                              min_types = 2) {
  # unify proximity windows -> "proximity"
  summary_df2 <- summary_df
  summary_df2[[type_col]][summary_df2[[type_col]] %in% c("win5", "win8", "win15")] <- "proximity"
  
  up <- summary_df2 %>%
    filter(.data[[type_col]] %in% c("colocal", "adjacency", "proximity")) %>%
    filter(.data[[score_col]] >= score_thresh, .data[[prop_col]] >= prop_thresh) %>%
    group_by(.data[[pair_col]]) %>%
    summarise(n_types_passing = n_distinct(.data[[type_col]]), .groups = "drop") %>%
    filter(n_types_passing >= min_types) %>%
    pull(.data[[pair_col]])
  
  down <- summary_df2 %>%
    filter(.data[[type_col]] %in% c("colocal", "adjacency", "proximity")) %>%
    filter(.data[[score_col]] <= -score_thresh, .data[[prop_col]] >= prop_thresh) %>%
    group_by(.data[[pair_col]]) %>%
    summarise(n_types_passing = n_distinct(.data[[type_col]]), .groups = "drop") %>%
    filter(n_types_passing >= min_types) %>%
    pull(.data[[pair_col]])
  
  list(up = up, down = down)
}

# Read the summary tables and format
# Output columns:
#   rownames, type (win5/win8/win15), prop_dep, prop_enr, pair_mean_rescaled, analysis
read_summary <- function(path, analysis_label) {
  x <- readRDS(path)
  
  if ("pair" %in% colnames(x) && !"rownames" %in% colnames(x)) {
    x <- dplyr::rename(x, rownames = pair)
  }
  
  # Expected new format (from your pipeline):
  # rownames, analysis (win5/win8/win15), prop_dep, prop_enr, mean_scaled
  if (all(c("rownames","analysis","prop_dep","prop_enr","mean_scaled") %in% colnames(x))) {
    x <- x %>%
      transmute(
        rownames = rownames,
        type = analysis,
        prop_dep = prop_dep,
        prop_enr = prop_enr,
        pair_mean_rescaled = mean_scaled
      )
  } else if (all(c("rownames","type","prop_dep","prop_enr","pair_mean_rescaled") %in% colnames(x))) {
    x <- x[, c("rownames","type","prop_dep","prop_enr","pair_mean_rescaled")]
  } else {
    stop("Unexpected summary table format in: ", path)
  }
  
  x$analysis <- analysis_label
  x$rownames <- gsub("\\bcHypResp\\b", "cMESVasc", x$rownames)
  x$rownames <- sort_pair(x$rownames)
  x
}


# Load colocal + adjacency and combine as "st_coloc_adj"


# ---- colocal ----
colocal_dir <- "/colocal/idh_gbm_comb_v2"

colocal_files <- c(
  low   = "colocal_catgrade_low.rds",
  med   = "colocal_catgrade_mid.rds",
  high  = "colocal_catgrade_high.rds",
  oligo = "colocal_oligo.rds",
  astro = "colocal_astro.rds",
  gbm   = "colocal_gbm.rds"
)

colocal_list <- imap(colocal_files, ~{
  df <- readRDS(file.path(colocal_dir, .x))
  df$analysis <- .y
  df
})

colocal_summary <- bind_rows(colocal_list)

colocal_summary$sig <- ifelse(
  colocal_summary$prop >= 0.25 & colocal_summary$pair_mean > 0, "yes_up",
  ifelse(colocal_summary$prop >= 0.25 & colocal_summary$pair_mean < 0, "yes_down", "no")
)

# rescale within analysis using pair_mean_scaled
colocal_summary <- rescale_signed_by_group(colocal_summary, value_col = "pair_mean_scaled", out_col = "pair_mean_rescaled")
colocal_summary <- colocal_summary[, c("rownames","analysis","prop","pair_mean","sig","type","pair_mean_rescaled")]

# make pair names consistent
colocal_summary$rownames <- str_replace_all(colocal_summary$rownames, "\\bcHypResp\\b", "cMESVasc")
colocal_summary$rownames <- sort_pair(colocal_summary$rownames)

# ---- adjacency ----
adj_dir <- "/adjacency/idh_gbm_comb_v2"

adj_files <- c(
  astro = "sum_adj_fin_astro.rds",
  oligo = "sum_adj_fin_oligo.rds",
  gbm   = "sum_adj_fin_catgrade_gbm.rds",
  low   = "sum_adj_fin_catgrade_low.rds",
  med   = "sum_adj_fin_catgrade_med.rds",
  high  = "sum_adj_fin_catgrade_high.rds"
)

adj_list <- imap(adj_files, ~{
  df <- readRDS(file.path(adj_dir, .x))
  df$analysis <- .y
  df
})

adj_summary <- bind_rows(adj_list)

adj_summary$rownames <- rownames(adj_summary)
adj_summary$type <- "adjacency"
adj_summary$sig <- ifelse(
  adj_summary$prop >= 0.25 & adj_summary$mean_score > 0, "yes_up",
  ifelse(adj_summary$prop >= 0.25 & adj_summary$mean_score < 0, "yes_down", "no")
)
adj_summary$pair_mean <- adj_summary$mean_score

adj_summary <- adj_summary[, c("rownames","analysis","prop","pair_mean","sig","type")]

# rescale 
adj_summary <- rescale_signed_by_group(adj_summary, value_col = "pair_mean", out_col = "pair_mean_rescaled")

# make pair names consistent
adj_summary$rownames <- str_replace_all(adj_summary$rownames, "\\bcHypResp\\b", "cMESVasc")
adj_summary$rownames <- sort_pair(adj_summary$rownames)

# combined coloc + adj
st_coloc_adj <- bind_rows(
  colocal_summary[, c("rownames","analysis","prop","pair_mean","sig","type","pair_mean_rescaled")],
  adj_summary[,   c("rownames","analysis","prop","pair_mean","sig","type","pair_mean_rescaled")]
)


#Load regional-comp summary tables and map into "summary_*" dfs


new_sum_dir <- "/reg_comp/summary_tables"

# TYPE summaries
summary_astro <- read_summary(file.path(new_sum_dir, "summary_type_astro.rds"), analysis_label = "astro")
summary_oligo <- read_summary(file.path(new_sum_dir, "summary_type_oligo.rds"), analysis_label = "oligo")
summary_gbm   <- read_summary(file.path(new_sum_dir, "summary_type_gbm.rds"),   analysis_label = "gbm")

# CategoricalGrade summaries
summary_low  <- read_summary(file.path(new_sum_dir, "summary_catgrade_low.rds"),  analysis_label = "low")
summary_med  <- read_summary(file.path(new_sum_dir, "summary_catgrade_med.rds"),  analysis_label = "med")
summary_high <- read_summary(file.path(new_sum_dir, "summary_catgrade_high.rds"), analysis_label = "high")

# Standardize regional-comp dfs to match downstream use
process_reg_df <- function(df) {
  df$sig <- ifelse(df$prop_enr >= 0.25, "yes_up",
                   ifelse(df$prop_dep >= 0.25, "yes_down", "no"))
  df$pair_mean <- df$pair_mean_rescaled
  df <- df[, c("rownames","analysis","prop_enr","prop_dep","pair_mean","sig","type","pair_mean_rescaled")]
  df
}

summary_astro <- process_reg_df(summary_astro)
summary_oligo <- process_reg_df(summary_oligo)
summary_gbm   <- process_reg_df(summary_gbm)
summary_low   <- process_reg_df(summary_low)
summary_med   <- process_reg_df(summary_med)
summary_high  <- process_reg_df(summary_high)


# Build per-group "summary_df" and compute consensus up/down


make_group_summary_df <- function(group,
                                  st_coloc_adj,
                                  reg_summary,
                                  pair_col = "rownames") {
  # colocal rows
  coloc_df <- st_coloc_adj %>%
    filter(analysis == group, type == "colocal") %>%
    transmute(rownames = .data[[pair_col]],
              type = "colocal",
              prop = prop,
              mean_scaled = pair_mean_rescaled)
  
  # adjacency rows
  adj_df <- st_coloc_adj %>%
    filter(analysis == group, type == "adjacency", sig != "no") %>%
    transmute(rownames = .data[[pair_col]],
              type = "adjacency",
              prop = prop,
              mean_scaled = pair_mean_rescaled)
  
  # proximity rows (regional comp)
  prox_df <- reg_summary %>%
    mutate(prop = prop_dep + prop_enr) %>%
    transmute(rownames = .data[[pair_col]],
              type = type,                 # win5/win8/win15
              prop = prop,
              mean_scaled = pair_mean_rescaled)
  
  bind_rows(coloc_df, adj_df, prox_df)
}

# TYPE groups
ast_summary_df   <- make_group_summary_df("astro", st_coloc_adj, summary_astro)
oligo_summary_df <- make_group_summary_df("oligo", st_coloc_adj, summary_oligo)
gbm_summary_df   <- make_group_summary_df("gbm",   st_coloc_adj, summary_gbm)

# GRADE groups
low_summary_df   <- make_group_summary_df("low",   st_coloc_adj, summary_low)
med_summary_df   <- make_group_summary_df("med",   st_coloc_adj, summary_med)
high_summary_df  <- make_group_summary_df("high",  st_coloc_adj, summary_high)

# compute consensus
cons_ast  <- compute_consensus(ast_summary_df,  pair_col = "rownames")
cons_olig <- compute_consensus(oligo_summary_df, pair_col = "rownames")
cons_gbm  <- compute_consensus(gbm_summary_df,  pair_col = "rownames")

cons_low  <- compute_consensus(low_summary_df,  pair_col = "rownames")
cons_med  <- compute_consensus(med_summary_df,  pair_col = "rownames")
cons_high <- compute_consensus(high_summary_df, pair_col = "rownames")

consensus_up_ast    <- cons_ast$up
consensus_down_ast  <- cons_ast$down

consensus_up_oligo   <- cons_olig$up
consensus_down_oligo <- cons_olig$down

consensus_up_gbm    <- cons_gbm$up
consensus_down_gbm  <- cons_gbm$down

consensus_up_low    <- cons_low$up
consensus_down_low  <- cons_low$down

consensus_up_med    <- cons_med$up
consensus_down_med  <- cons_med$down

consensus_up_high   <- cons_high$up
consensus_down_high <- cons_high$down


#consensus interactions

cons <- list(
  "consensus_gbm"   = consensus_up_gbm,
  "consensus_astro" = consensus_up_ast,
  "consensus_oligo" = consensus_up_oligo,
  "consensus_low"   = consensus_up_low,
  "consensus_med"   = consensus_up_med,
  "consensus_high"  = consensus_up_high,
  "repel_gbm"       = consensus_down_gbm,
  "repel_astro"     = consensus_down_ast,
  "repel_oligo"     = consensus_down_oligo,
  "repel_low"       = consensus_down_low,
  "repel_med"       = consensus_down_med,
  "repel_high"      = consensus_down_high
)

saveRDS(cons, "consensus_interactions_list.rds")


