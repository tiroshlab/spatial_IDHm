library(dplyr)
library(ggplot2)
library(patchwork)


##########Interaction type analysis####################

cons<-readRDS("consensus_interactions_list.rds")

######Defining consensus pairs and interaction types -- Allowed pairs across all lists ----------


allowed_pairs <- Reduce(union, list(
  cons$consensus_gbm, cons$consensus_astro, cons$consensus_oligo,
  cons$repel_gbm,     cons$repel_astro,     cons$repel_oligo
))

interaction_type <- unname(
  sapply(allowed_pairs, function(p) {
    states <- strsplit(p, " ")[[1]]
    paste(ifelse(grepl("^c", states), "Cancer", "TME"), collapse = "-")
  })
)


interaction_type_df<-cbind(allowed_pairs,interaction_type)
interaction_type_df->allowed_pairs_df

rename_interaction <- function(x) dplyr::recode(
  x,
  "malignant_malignant"        = "Cancer-Cancer",
  "malignant_nonmalignant"     = "Cancer-TME",
  "nonmalignant_nonmalignant"  = "TME-TME",
  .default = x
)

###### Interaction heatmap by type


# Order columns (allowed_pairs) by fixed interaction-type blocks
#    in the order: TME-TME -> Cancer-TME -> Cancer-Cancer,
#    and within each block, put GBM-present pairs first, then alphabetically.
row_order_by_fixed_type_gbm_first <- function(pairs_union, gbm_pairs, lookup_df) {
  desired <- c("TME-TME", "Cancer-TME", "Cancer-Cancer")
  
  sub_df <- lookup_df %>%
    dplyr::filter(allowed_pairs %in% pairs_union) %>%
    dplyr::mutate(
      interaction_type = factor(rename_interaction(interaction_type), levels = desired),
      in_gbm = allowed_pairs %in% gbm_pairs
    )
  
  sub_df %>%
    dplyr::arrange(interaction_type, dplyr::desc(in_gbm), allowed_pairs) %>%
    dplyr::pull(allowed_pairs) %>%
    unique()
}


make_binary_long <- function(row_order, pairs_list, labels, lookup_df) {
  ref_grid <- expand.grid(
    allowed_pairs = row_order,
    cohort = labels,
    stringsAsFactors = FALSE
  )
  
  present_df <- do.call(rbind, lapply(seq_along(pairs_list), function(i) {
    data.frame(
      allowed_pairs = unique(pairs_list[[i]]),
      cohort = labels[i],
      present = TRUE,
      stringsAsFactors = FALSE
    )
  }))
  
  ref_grid %>%
    dplyr::left_join(present_df, by = c("allowed_pairs","cohort")) %>%
    dplyr::mutate(present = ifelse(is.na(present), FALSE, present)) %>%
    dplyr::left_join(lookup_df, by = "allowed_pairs") %>%
    dplyr::mutate(
      fill_type     = ifelse(present, interaction_type, NA_character_),
      allowed_pairs = factor(allowed_pairs, levels = row_order),
      cohort        = factor(cohort, levels = labels)
    )
}

plot_binary <- function(df_long, title_txt, fill_vals) {
  df_long <- df_long %>%
    dplyr::mutate(fill_type = dplyr::recode(
      fill_type,
      "malignant_malignant"        = "Cancer-Cancer",
      "malignant_nonmalignant"     = "Cancer-TME",
      "nonmalignant_nonmalignant"  = "TME-TME",
      .default = fill_type
    ))
  
  ggplot(df_long, aes(y = cohort, x = allowed_pairs, fill = fill_type)) +
    geom_tile(color = "grey90", linewidth = 0.4) +
    scale_fill_manual(values = fill_vals, na.value = "white", drop = FALSE,
                      name = "Interaction category") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 0, vjust = 0, size = 10),
      axis.text.y  = element_text(size = 14),
      plot.title   = element_text(hjust = 0.5, size = 18, margin = margin(b = 10)),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.margin  = margin(t = 30, r = 10, b = 10, l = 10),
      legend.position = "right"
    )
}

pal <- c(
  "Cancer-Cancer" = "#ff0066",
  "Cancer-TME"    = "#a5506d",
  "TME-TME"       = "#328c97"
)

# Ensure lookup is a data.frame (if it isn't already)
allowed_pairs_df <- as.data.frame(allowed_pairs_df)

# --- CONSENSUS panel ---
consensus_union <- Reduce(union, list(cons$consensus_gbm, cons$consensus_astro, cons$consensus_oligo))
consensus_union <- intersect(consensus_union, allowed_pairs_df$allowed_pairs)

row_order_consensus <- row_order_by_fixed_type_gbm_first(
  pairs_union = consensus_union,
  gbm_pairs   = intersect(consensus_union, cons$consensus_gbm),
  lookup_df   = allowed_pairs_df
)

consensus_long <- make_binary_long(
  row_order   = row_order_consensus,
  pairs_list  = list(cons$consensus_gbm, cons$consensus_astro, cons$consensus_oligo),
  labels      = c("GBM","IDH-A","IDH-O"),
  lookup_df   = allowed_pairs_df
)

p_consensus <- plot_binary(consensus_long, "Consensus interaction", pal)



# stacked barplots - Fig. 3D interaction type proportions by tumor type -----------

########proportion of interaction type per cancer type - consensus + repelled
rename_interaction <- function(x) {
  dplyr::recode(
    x,
    "malignant_malignant"      = "Cancer-Cancer",
    "malignant_nonmalignant"   = "Cancer-TME",
    "nonmalignant_nonmalignant"= "TME-TME",
    .default = x
  )
}

## helper: interaction-type order (GBM-first), with renamed labels
type_order_by_gbm <- function(pairs_union, gbm_pairs, lookup_df) {
  sub_df <- lookup_df %>% dplyr::filter(allowed_pairs %in% pairs_union)
  
  gbm_types <- sub_df %>%
    dplyr::filter(allowed_pairs %in% gbm_pairs) %>%
    dplyr::pull(interaction_type) %>%
    rename_interaction() %>%
    unique()
  
  remaining_types <- setdiff(
    sort(unique(rename_interaction(sub_df$interaction_type))),
    gbm_types
  )
  
  c(gbm_types, remaining_types)
}

## build stacked bars (horizontal), using global `pal`
make_stacked_bar <- function(df_long, pairs_union, gbm_pairs, lookup_df, title_txt) {
  # GBM-first stack order (renamed)
  type_order <- type_order_by_gbm(pairs_union, gbm_pairs, lookup_df)
  
  bar_df <- df_long %>%
    dplyr::filter(present, !is.na(fill_type)) %>%
    dplyr::mutate(interaction_type = factor(rename_interaction(fill_type),
                                            levels = type_order)) %>%
    dplyr::count(cohort, interaction_type, name = "n") %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  ggplot(bar_df, aes(x = cohort, y = n, fill = interaction_type)) +
    geom_bar(stat = "identity", position = "fill") +     # proportions per cohort
    coord_flip() +                                        # <-- horizontal bars
    scale_fill_manual(values = pal, drop = FALSE, name = "Interaction type") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion", title = title_txt) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x  = element_text(size = 14, angle=90, hjust=1, vjust=0.5),
      axis.text.y  = element_text(size = 14),
      plot.title   = element_text(hjust = 0.5, size = 18, margin = margin(b = 8)),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      legend.position = "right"
    )
}

# calls (unchanged inputs; no fill_vals arg needed)

type_order_consensus <- type_order_by_gbm(
  pairs_union = consensus_union,
  gbm_pairs   = intersect(consensus_union, cons$consensus_gbm),
  lookup_df   = allowed_pairs_df
)
p_consensus_bar <- make_stacked_bar(
  df_long     = consensus_long,
  pairs_union = consensus_union,
  gbm_pairs   = intersect(consensus_union, cons$consensus_gbm),
  lookup_df   = allowed_pairs_df,
  title_txt   = "Proportion interaction category"
)



# line plot - consensus int and repelled int count by tumor type: Fig 3B

#consensus
consensus_counts <- data.frame(
  CancerType = c("GBM", "IDH-A", "IDH-O"),
  Count = c(length(cons$consensus_gbm),
            length(cons$consensus_astro),
            length(cons$consensus_oligo))
)

# simple line plot
p1 <- ggplot(consensus_counts, aes(x = CancerType, y = Count, group = 1)) +
  geom_line(color = "#B01326", size = 1.5) +
  geom_point(color = "#B01326", size = 4) +
  labs(x = "Glioma type", y = "# consensus interactions") +
  theme_classic(base_size = 14)+
  expand_limits(y = 0)

#repelled
repelled_counts <- data.frame(
  CancerType = c("GBM", "IDH-A", "IDH-O"),
  Count = c(length(cons$repel_gbm),
            length(cons$repel_astro),
            length(cons$repel_oligo))
)

# simple line plot
p2 <- ggplot(repelled_counts, aes(x = CancerType, y = Count, group = 1)) +
  geom_line(color = "#2166AC", size = 1.5) +
  geom_point(color = "#2166AC", size = 4) +
  labs(x = "Glioma type", y = "# consensus repulsions") +
  theme_classic(base_size = 14)+
  expand_limits(y = 0)

#prop. consensus vs. repelled by type

counts <- data.frame(
  CancerType = rep(c("GBM","IDH-A","IDH-O"), each = 2),
  Category   = rep(c("Consensus","Repelled"), times = 3),
  Count = c(length(cons$consensus_gbm), length(cons$repel_gbm),
            length(cons$consensus_astro), length(cons$repel_astro),
            length(cons$consensus_oligo), length(cons$repel_oligo))
)

p1/p2




#########interaction type plots by progression level##############

####Allowed pairs across all lists ----------
cons$consensus_low->consensus_low
cons$consensus_med->consensus_med
cons$consensus_high->consensus_high
cons$consensus_gbm->consensus_gbm

cons$repel_low->repel_low
cons$repel_med->repel_med
cons$repel_high->repel_high
cons$repel_gbm->repel_gbm


allowed_pairs <- Reduce(union, list(
  consensus_low, consensus_med, consensus_high, consensus_gbm,
  repel_low, repel_med, repel_high, repel_gbm
))

interaction_type <- unname(
  sapply(allowed_pairs, function(p) {
    states <- strsplit(p, " ")[[1]]
    paste(ifelse(grepl("^c", states), "Cancer", "TME"), collapse = "-")
  })
)

interaction_type_df<-cbind(allowed_pairs,interaction_type)


#######interaction type binary heatmaps - by grade

interaction_type_df <- as.data.frame(interaction_type_df, stringsAsFactors = FALSE)


stopifnot(all(c("allowed_pairs","interaction_type") %in% colnames(interaction_type_df)))
interaction_type_df$allowed_pairs    <- as.character(interaction_type_df$allowed_pairs)
interaction_type_df$interaction_type <- as.character(interaction_type_df$interaction_type)


allowed_pairs_df <- interaction_type_df

## ---- helper functions
row_order_by_low <- function(pairs_union, low_pairs, lookup_df) {
  # Ensure proper type
  lookup_df <- as.data.frame(lookup_df, stringsAsFactors = FALSE)
  if ("allowed_pairs" %in% names(lookup_df))    lookup_df$allowed_pairs    <- as.character(lookup_df$allowed_pairs)
  if ("interaction_type" %in% names(lookup_df)) lookup_df$interaction_type <- as.character(lookup_df$interaction_type)
  
  sub_df <- lookup_df %>% dplyr::filter(allowed_pairs %in% pairs_union)
  
  # order interaction types by the order they appear among 'low' pairs (first-seen),
  # then append any remaining types alphabetically
  low_types <- sub_df %>%
    dplyr::filter(allowed_pairs %in% low_pairs) %>%
    dplyr::pull(interaction_type) %>%
    unique()
  
  remaining_types <- setdiff(sort(unique(sub_df$interaction_type)), low_types)
  type_order <- c(low_types, remaining_types)
  
  sub_df %>%
    dplyr::mutate(interaction_type = factor(interaction_type, levels = type_order)) %>%
    dplyr::arrange(interaction_type, allowed_pairs) %>%
    dplyr::pull(allowed_pairs) %>%
    unique()
}

# Build long binary df for N grade columns with presence colored by interaction_type
make_binary_long <- function(row_order, pairs_list, labels, lookup_df) {
  # Ensure proper type
  lookup_df <- as.data.frame(lookup_df, stringsAsFactors = FALSE)
  if ("allowed_pairs" %in% names(lookup_df))    lookup_df$allowed_pairs    <- as.character(lookup_df$allowed_pairs)
  if ("interaction_type" %in% names(lookup_df)) lookup_df$interaction_type <- as.character(lookup_df$interaction_type)
  
  ref_grid <- expand.grid(
    allowed_pairs = row_order,
    cohort = labels,
    stringsAsFactors = FALSE
  )
  
  present_df <- do.call(rbind, lapply(seq_along(pairs_list), function(i) {
    data.frame(
      allowed_pairs = unique(pairs_list[[i]]),
      cohort = labels[i],
      present = TRUE,
      stringsAsFactors = FALSE
    )
  }))
  
  dplyr::left_join(ref_grid, present_df, by = c("allowed_pairs","cohort")) %>%
    dplyr::mutate(present = ifelse(is.na(present), FALSE, present)) %>%
    dplyr::left_join(lookup_df, by = "allowed_pairs") %>%
    dplyr::mutate(
      fill_type    = ifelse(present, interaction_type, NA_character_),
      allowed_pairs = factor(allowed_pairs, levels = row_order),
      cohort        = factor(cohort, levels = labels)
    )
}

#####helper function for interaction type binary heatmaps by interaction type
plot_binary <- function(df_long, title_txt, fill_vals) {
  df_long <- df_long %>%
    mutate(fill_type = recode(
      fill_type,
      "malignant_malignant" = "Cancer-Cancer",
      "malignant_nonmalignant" = "Cancer-TME",
      "nonmalignant_nonmalignant" = "TME-TME"
    ))
  
  ggplot(df_long, aes(y = cohort, x = allowed_pairs, fill = fill_type)) +
    geom_tile(color = "grey90", linewidth = 0.4) +
    scale_fill_manual(values = fill_vals, na.value = "white", drop = FALSE,
                      name = "Interaction type") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 0, vjust = 0, size = 14),
      axis.text.y   = element_text(size = 14),
      plot.title    = element_text(hjust = 0.5, size = 18, margin = margin(b = 10)),
      legend.title  = element_text(size = 16),
      legend.text   = element_text(size = 14),
      plot.margin   = margin(t = 30, r = 10, b = 10, l = 10),
      legend.position = "right"
    )
}


interaction_types_all <- sort(unique(allowed_pairs_df$interaction_type))

pal <- c(
  "Cancer-Cancer" = "#ff0066",
  "Cancer-TME" =  "#a5506d",
  "TME-TME" = "#328c97"
)

## labels 
grade_labels <- c("Low","Mid","High","GBM")

#######interaction type binary heatmaps by interaction type - Fig. 3G#########

## CONSENSUS panel 
consensus_union <- Reduce(union, list(consensus_low, consensus_med, consensus_high, consensus_gbm))
consensus_union <- intersect(consensus_union, allowed_pairs_df$allowed_pairs)

row_order_consensus <- row_order_by_low(
  pairs_union = consensus_union,
  low_pairs   = intersect(consensus_union, consensus_low),
  lookup_df   = allowed_pairs_df
)

consensus_long <- make_binary_long(
  row_order  = row_order_consensus,
  pairs_list = list(consensus_low, consensus_med, consensus_high, consensus_gbm),
  labels     = grade_labels,
  lookup_df  = allowed_pairs_df
)

p_consensus <- plot_binary(consensus_long, "Consensus interaction", pal)

## REPELLED panel
repelled_union <- Reduce(union, list(repel_low, repel_med, repel_high, repel_gbm))
repelled_union <- intersect(repelled_union, allowed_pairs_df$allowed_pairs)

row_order_repelled <- row_order_by_low(
  pairs_union = repelled_union,
  low_pairs   = intersect(repelled_union, repel_low),
  lookup_df   = allowed_pairs_df
)

repelled_long <- make_binary_long(
  row_order  = row_order_repelled,
  pairs_list = list(repel_low, repel_med, repel_high, repel_gbm),
  labels     = grade_labels,
  lookup_df  = allowed_pairs_df
)

p_repelled <- plot_binary(repelled_long, "Consensus repulsion", pal)

p_consensus / p_repelled +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

########proportion of interaction type  - stacked barplot consensus + repelled - fig 3H
##interaction-type order 
type_order_by_low <- function(pairs_union, low_pairs, lookup_df) {
  sub_df <- lookup_df %>% dplyr::filter(allowed_pairs %in% pairs_union)
  low_types <- sub_df %>%
    dplyr::filter(allowed_pairs %in% low_pairs) %>%
    dplyr::pull(interaction_type) %>%
    unique()
  remaining_types <- setdiff(sort(unique(sub_df$interaction_type)), low_types)
  c(low_types, remaining_types)
}

## ---- build stacked bars --------------------------------------------------
make_stacked_bar <- function(df_long, pairs_union, low_pairs, lookup_df, fill_vals, title_txt) {
  # GBM-first stack order
  type_order <- type_order_by_low(pairs_union, low_pairs, lookup_df)
  
  bar_df <- df_long %>%
    dplyr::filter(present, !is.na(fill_type)) %>%
    dplyr::mutate(interaction_type = factor(fill_type, levels = type_order)) %>%
    dplyr::count(cohort, interaction_type, name = "n") %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  ggplot(bar_df, aes(x = cohort, y = n, fill = interaction_type)) +
    geom_bar(stat = "identity", position = "fill") +  # proportions per cohort
    scale_fill_manual(
      values = pal,
      drop = FALSE,
      name = "Interaction type"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion", title = title_txt) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
      axis.text.y = element_text(size = 16),
      plot.title  = element_text(hjust = 0.5, size = 18, margin = margin(b = 8)),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 16),
      legend.position = "none"
    )
  
}

## ---- consensus stacked bar ----------------------------------------------
type_order_consensus <- type_order_by_low(
  pairs_union = consensus_union,
  low_pairs   = intersect(consensus_union, consensus_low),
  lookup_df   = allowed_pairs_df
)
p_consensus_bar <- make_stacked_bar(
  df_long     = consensus_long,
  pairs_union = consensus_union,
  low_pairs   = intersect(consensus_union, consensus_low),
  lookup_df   = allowed_pairs_df,
  fill_vals   = pal,
  title_txt   = "Consensus interactions: proportion by grade"
)

## ---- repelled stacked bar -----------------------------------------------
type_order_repelled <- type_order_by_low(
  pairs_union = repelled_union,
  low_pairs   = intersect(repelled_union, repel_low),
  lookup_df   = allowed_pairs_df
)
p_repelled_bar <- make_stacked_bar(
  df_long     = repelled_long,
  pairs_union = repelled_union,
  low_pairs   = intersect(repelled_union, repel_low),
  lookup_df   = allowed_pairs_df,
  fill_vals   = fill_vals,
  title_txt   = "Repelled interactions: proportion by grade"
)

## ---- show ---------------------------------------------------------------
p_consensus_bar + p_repelled_bar


####line plot - consensus int and repelled int count by prog level - fig 3f 

#consensus
consensus_counts <- data.frame(
  CancerType = c("low-grade", "med-grade", "high-grade","GBM"),
  Count = c(length(consensus_low),
            length(consensus_med),
            length(consensus_high),
            length(consensus_gbm))
)

# set factor order
consensus_counts$CancerType <- factor(
  consensus_counts$CancerType,
  levels = c("low-grade", "med-grade", "high-grade","GBM")
)


#repelled
repelled_counts <- data.frame(
  CancerType = c("low-grade", "med-grade", "high-grade","GBM"),
  Count = c(length(repel_low),
            length(repel_med),
            length(repel_high),
            length(repel_gbm))
)

# set factor order
repelled_counts$CancerType <- factor(
  repelled_counts$CancerType,
  levels = c("low-grade", "med-grade", "high-grade","GBM"))


p1 <- ggplot(consensus_counts, aes(x = CancerType, y = Count, group = 1)) +
  geom_line(color = "#B01326", size = 1.5) +
  geom_point(color = "#B01326", size = 4) +
  labs(x = "Glioma type", y = "# consensus interactions") +
  theme_classic(base_size = 14) +
  expand_limits(y = 0)

p2 <- ggplot(repelled_counts, aes(x = CancerType, y = Count, group = 1)) +
  geom_line(color = "#2166AC", size = 1.5) +
  geom_point(color = "#2166AC", size = 4) +
  labs(x = "Glioma type", y = "# consensus repulsions") +
  theme_classic(base_size = 14) +
  expand_limits(y = 0)

p1/p2

###################prop. consensus vs. repelled by type - Fig S4B/C

# Stacked bars: proportion of relationship type (interaction vs repulsion)
#   - by glioma type (GBM / IDH-A / IDH-O)
#   - by progression/grade (Low / Mid / High / GBM)


library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

# Colors to match your example
rel_pal <- c(
  "interaction" = "#B01326",  # red
  "repulsion"   = "#2166AC"   # blue
)

make_rel_prop_plot <- function(df_counts, xvar, title_txt = NULL) {
  df_long <- df_counts %>%
    tidyr::pivot_longer(c(interaction, repulsion),
                        names_to = "relationship", values_to = "n") %>%
    group_by(.data[[xvar]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(relationship = factor(relationship, levels = c("repulsion","interaction"))) # blue bottom
  
  ggplot(df_long, aes(x = .data[[xvar]], y = prop, fill = relationship)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = rel_pal, name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion of\nrelationship type", title = title_txt) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = "right",
      legend.text = element_text(size = 14)
    )
}


# By glioma TYPE (GBM / IDH-A / IDH-O)

counts_type <- tibble(
  GliomaType  = factor(c("GBM","IDH-A","IDH-O"), levels = c("GBM","IDH-A","IDH-O")),
  interaction = c(length(cons$consensus_gbm),  length(cons$consensus_astro), length(cons$consensus_oligo)),
  repulsion   = c(length(cons$repel_gbm),      length(cons$repel_astro),     length(cons$repel_oligo))
)

p_rel_type <- make_rel_prop_plot(counts_type, "GliomaType", "Glioma type")


#By progression / GRADE (Low / Mid / High / GBM)

counts_grade <- tibble(
  Grade       = factor(c("Low","Mid","High","GBM"), levels = c("Low","Mid","High","GBM")),
  interaction = c(length(consensus_low), length(consensus_med), length(consensus_high), length(consensus_gbm)),
  repulsion   = c(length(repel_low),     length(repel_med),     length(repel_high),     length(repel_gbm))
)

p_rel_grade <- make_rel_prop_plot(counts_grade, "Grade", "Progression / grade")

# Show (side-by-side or stacked)
p_rel_type | p_rel_grade
# or:
# p_rel_type / p_rel_grade


