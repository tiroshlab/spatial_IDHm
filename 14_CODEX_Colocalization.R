# Custom functions --------------------------------------------------------
marker <- rownames(rmat)

## Test neighborhood with shuffling
Test_nhood <- function(sample = "IDHA01", 
                       radius = 27.5, 
                       cell_type_column = "cell_type", 
                       fraction_coherance = 0.8,
                       iter = 10,
                       workers = 10) {
  
  # Observed part  
  cells_clean <- cells %>% 
    filter(cell_type %ni% c("excluded", "low")) %>% 
    filter(sample %in% c({{sample}}))
  
  spe_clean <- SpatialExperiment::SpatialExperiment(
    sample_id = cells_clean$sample,
    spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
    rowData = marker,
    colData = list(cell_name = cells_clean$cell_name,
                   cell_type = as.factor(cells_clean[[cell_type_column]]),
                   centroid_x = cells_clean$centroid_x,
                   centroid_y = cells_clean$centroid_y,
                   img_id = cells_clean$sample,
                   all = rep("all", dim(cells_clean)[1]))
  )
  
  spe_clean <- imcRtools::buildSpatialGraph(spe_clean, 
                                            img_id = "sample_id", 
                                            type = "expansion", 
                                            threshold = radius, 
                                            coords = c("centroid_x", "centroid_y"),
                                            name = "expansion_graph")
  
  
  spe_clean <- imcRtools::aggregateNeighbors(spe_clean, 
                                             colPairName = "expansion_graph", 
                                             aggregate_by = "metadata", 
                                             count_by = "cell_type",
                                             name = "nhood_mat", 
                                             proportions = F)
  
  
  nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                     from_cell_type = cells_clean[[cell_type_column]]) %>% as_tibble
  
  
  nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
  
  row_sums <- nhood_mat_num %>% rowSums()
  fraction_coherance_mat <- nhood_mat_num/row_sums
  drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance))
  nhood_norm_coherant <- nhood_mat[!drop_rows,]
  nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
  
  ## Summarise interactions
  
  nhood_sum <- nhood_norm_coherant %>% 
    group_by(from_cell_type) %>% 
    summarise(across(where(is.numeric), sum))
  
  ## Calc rowSums and subtract self-pairs
  nhood_rowSums <- nhood_sum %>% 
    select_if(is.numeric) %>% 
    rowSums()
  
  ## Pivot_long and extract self-pairs
  self_pair_count <- nhood_sum %>% 
    pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = "count") %>% 
    filter(from_cell_type == to_cell_type) %>% 
    pull(count)
  
  ## Normalize nhood_sum by difference rowSums and self_pair_count
  
  nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
    cbind(from_cell_type = nhood_sum$from_cell_type, .) %>% as_tibble()
  
  nhood_norm_obs <- nhood_norm %>% 
    pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = paste0("obs_count_" ,{{sample}}))
  
  # Shuffeled part 
  
  Tmp <- function(sample. = sample, 
                  radius. = radius, 
                  cell_type_column. = cell_type_column, 
                  fraction_coherance. = fraction_coherance,
                  iter. = iter,
                  workers. = workers) 
  {
    
    shuffled_counts <- BiocParallel::bplapply(1:iter., 
                                              BPPARAM = BiocParallel::MulticoreParam(workers = workers., progressbar = T),
                                              function(i) {
                                                
                                                cells_clean <- cells %>% 
                                                  filter(cell_type %ni% c("excluded", "low")) %>% 
                                                  filter(sample %in% c(sample.))
                                                
                                                cells_clean[[cell_type_column.]] <- sample(cells_clean[[cell_type_column.]])
                                                
                                                spe_clean <- SpatialExperiment::SpatialExperiment(
                                                  sample_id = cells_clean$sample,
                                                  spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
                                                  rowData = marker,
                                                  colData = list(cell_name = cells_clean$cell_name,
                                                                 cell_type = as.factor(cells_clean[[cell_type_column.]]),
                                                                 centroid_x = cells_clean$centroid_x,
                                                                 centroid_y = cells_clean$centroid_y,
                                                                 img_id = cells_clean$sample,
                                                                 all = rep("all", dim(cells_clean)[1]))
                                                )
                                                
                                                spe_clean <- imcRtools::buildSpatialGraph(spe_clean, 
                                                                                          img_id = "sample_id", 
                                                                                          type = "expansion", 
                                                                                          threshold = radius., 
                                                                                          coords = c("centroid_x", "centroid_y"),
                                                                                          name = "expansion_graph")
                                                
                                                
                                                spe_clean <- imcRtools::aggregateNeighbors(spe_clean, 
                                                                                           colPairName = "expansion_graph", 
                                                                                           aggregate_by = "metadata", 
                                                                                           count_by = "cell_type",
                                                                                           name = "nhood_mat", 
                                                                                           proportions = F)
                                                
                                                
                                                nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                                                                   from_cell_type = cells_clean[[cell_type_column.]]) %>% as_tibble
                                                
                                                
                                                nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
                                                
                                                row_sums <- nhood_mat_num %>% rowSums()
                                                fraction_coherance_mat <- nhood_mat_num/row_sums
                                                drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance.))
                                                nhood_norm_coherant <- nhood_mat[!drop_rows,]
                                                nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
                                                
                                                ## Summarise interactions
                                                
                                                nhood_sum <- nhood_norm_coherant %>% 
                                                  group_by(from_cell_type) %>% 
                                                  summarise(across(where(is.numeric), sum))
                                                
                                                ## Calc rowSums and subtract self-pairs
                                                nhood_rowSums <- nhood_sum %>% 
                                                  select_if(is.numeric) %>% 
                                                  rowSums()
                                                
                                                ## Pivot_long and extract self-pairs
                                                self_pair_count <- nhood_sum %>% 
                                                  pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = "count") %>% 
                                                  filter(from_cell_type == to_cell_type) %>% 
                                                  pull(count)
                                                
                                                ## Normalize nhood_sum by difference rowSums and self_pair_count
                                                
                                                nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
                                                  cbind(from_cell_type = nhood_sum$from_cell_type, .) %>% as_tibble()
                                                
                                                nhood_norm_shuff <- nhood_norm %>% 
                                                  pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = paste0("shuff_", i)) #%>% 
                                                #select(paste0("shuff_", i))
                                                
                                                return(nhood_norm_shuff)
                                              }
    )
    
  }
  
  list_nhood_norm_shuff <- Tmp()
  
  nhood_norm_shuff <-
    left_join(nhood_norm_obs,
              reduce(
                list_nhood_norm_shuff,
                left_join,
                by = c("from_cell_type", "to_cell_type")
              )) %>% as_tibble()
  
  nhood_scaled <-
    nhood_norm_shuff %>% select_if(is.numeric) %>% t %>%  scale(center = T, scale = T) %>% t
  
  nhood_scaled <-
    cbind(
      nhood_norm_shuff %>% select(from_cell_type, to_cell_type),
      nhood_scaled
    ) %>% as_tibble
  
  
  return(nhood_scaled)
  
}

Summarise_samples_nhood <- function(df, p_thresh = 0.01) {
  df %>%
    # 1) rename the one obs_count_* column to a uniform name
    rename(obs_count = matches("^obs_count_")) %>%
    
    # 2) compute null‐distribution mean & sd, z_score, empirical p‐value, and interaction_type
    mutate(
      perm_mean = rowMeans(across(matches("^shuff_")), na.rm = TRUE),
      perm_sd   = apply(select(., matches("^shuff_")), 1, sd, na.rm = TRUE),
      z_score   = (obs_count - perm_mean) / perm_sd,
      p_value   = ( rowSums(
        abs(across(matches("^shuff_"))) >= abs(obs_count)
      ) + 1 ) / (ncol(select(., matches("^shuff_"))) + 1),
      interaction_type = if_else(z_score > 0, "attraction", "avoidance")
    ) %>%
    
    # 3) build the pair string and flag significance
    mutate(
      pair = paste0(from_cell_type, "_", to_cell_type),
      sig  = p_value < p_thresh
    ) %>%
    
    # 4) select only the columns you want
    select(
      from_cell_type,
      to_cell_type,
      z_score,
      p_value,
      sig,
      interaction_type,
      perm_mean,
      perm_sd,
      pair
    )
}

# Plotting inputs ----------------------------------------------------------------
list_nhood <- readRDS("inputs/colocalization_sample_list.RDS")
order_cell_type2 <- c(
  "Neuron",
  "Oligo",
  "Astro",
  "InfMg",
  "cOPC",
  "cAC1",
  "cAC2",
  "cUndiff",
  "cMES",
  "cMESHyp",
  "Neutro",
  "GAM",
  "MacScav",
  "Vasc",
  "MacBorder",
  "TcellCD4",
  "TcellCD8",
  "Bcell"
)

order_cell_type4 <- c(
  "Neuron",
  "Oligo",
  "cOPC",
  "cAC",
  "cUndiff",
  "cMES",
  "cMESHyp",
  "Mac",
  "Vasc",
  "excluded"
)


#plot 10*22
hotmap <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "white",
            "#f7f7f7", "white", "#f4a582", "#d6604d", "#b2182b", "#67001f")

# Corresponding values scaled to range [-15, 15] with white centered from -5 to +5
color_values <- scales::rescale(c(-15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15), to = c(0, 1))
# Co-localisation analysis ------------------------------------------------
## Run this code or alternatively import "colocalization_sample_list.RDS"
# cells <- cells %>%
#   mutate(cell_type_tmp = case_when(
#     cell_type == "Vasc" ~ "Vasc",
#     TRUE ~ cell_type2
#   ))
# 
# list_nhood <- lapply(cells$sample %>% unique(), function(x) Test_nhood(sample = x,
#                                                                        iter = 500,
#                                                                        workers = 20,
#                                                                        cell_type_column = "cell_type_tmp",
#                                                                        radius = 27.5,
#                                                                        fraction_coherance = 1))
# names(list_nhood) <- cells$sample %>% unique()

# Plot individual samples -------------------------------------------------
# Apply across your list:
summary_samples <- purrr::map(list_nhood, Summarise_samples_nhood)

summary_samples_flat <- imap_dfr(summary_samples, ~ mutate(.x, sample = .y))

## Remove pairs per sample for which one partner has < 20 cells in that sample
# 1) get cell counts per sample × cell_type_tmp
cell_counts <- cells %>%
  count(sample, cell_type4, name = "n_cells")

# 2) join onto your summary and drop low‐count rows
summary_samples_flat <- summary_samples_flat %>%
  # bring in n_from for from_cell_type
  left_join(cell_counts, 
            by = c("sample", "from_cell_type" = "cell_type4")) %>%
  rename(n_from = n_cells) %>%
  # bring in n_to for to_cell_type
  left_join(cell_counts, 
            by = c("sample", "to_cell_type"   = "cell_type4")) %>%
  rename(n_to = n_cells) %>%
  # replace missing with 0 (if a cell type never appears in that sample)
  mutate(
    n_from = coalesce(n_from, 0L),
    n_to   = coalesce(n_to,   0L)
  ) %>%
  # only keep rows where both > 20
  filter(n_from > 20, n_to > 20) %>%
  # optionally drop the helper cols
  select(-n_from, -n_to)

## Plot coloc per sample
#select samples
tmp <- summary_samples_flat %>% 
  filter(sample %in% samples_IDH_A)

ggplot(tmp %>%  
         filter(from_cell_type != to_cell_type) %>% 
         drop_na(),
       aes(x = factor(to_cell_type, levels = order_cell_type4), 
           y = factor(from_cell_type, levels = order_cell_type4), 
           col=z_score)) +
  # geom_tile(col = "black", size=0, fill = "white") +
  geom_point(size=4) + 
  scale_color_gradientn(colors = hotmap, values = color_values,
                        limits = c(-15, 15), oob = scales::squish, na.value = "white") +
  theme_classic() + 
  facet_wrap(~sample) +
  guides(size = guide_legend(override.aes = list(color = "black"))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1)) + 
  theme(axis.text.y =element_text()) +
  labs(x=NULL, y=NULL, size = "significant in \nfraction of samples") +
  ggtitle("IDH-A")

# Summarize cohort --------------------------------------------------------
## Calculate z_score summaries and asses fraction of samples where this pair is present and significant
# Step 1: Total number of samples per (from, to) pair
pair_totals <- summary_samples_flat %>%
  filter(sample %in% samples_IDH_A) %>% 
  distinct(sample, from_cell_type, to_cell_type) %>%
  count(from_cell_type, to_cell_type, name = "n_total")

# Step 2: Interaction-type-specific summary, with `n_significant`
summary_tbl <- summary_samples_flat %>%
  filter(sample %in% samples_IDH_A) %>% 
  group_by(from_cell_type, to_cell_type) %>%
  summarise(
    mean_z = mean(z_score, na.rm = TRUE),
    median_z = median(z_score, na.rm = TRUE),
    n_sig = sum(sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pair_totals, by = c("from_cell_type", "to_cell_type")) %>%
  mutate(percent_samples = n_sig / n_total)

## Plot
ggplot(summary_tbl %>%  
         filter(from_cell_type != to_cell_type) %>% 
         drop_na() %>% 
         filter(!from_cell_type %in% c("excluded"),
                !to_cell_type %in% c("excluded")),
       aes(x = factor(to_cell_type, levels = order_cell_type4), 
           y = factor(from_cell_type, levels = order_cell_type4), 
           col=median_z,
           size=percent_samples)) +
  geom_tile(col = "black", size=0, fill = "white") +
  geom_point(shape=19) + 
  scale_color_gradientn(colors = hotmap, values = color_values,
                        limits = c(-15, 15), oob = scales::squish, na.value = "white") +
  scale_size(range=c(5,12),limits = c(0.1,1)) +
  theme_classic() + 
  guides(size = guide_legend(override.aes = list(color = "black"))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16)) + 
  labs(x="Neighbor", y="Reference", size = "significant in \nfraction of samples") +
  ggtitle("IDH-A")

summary_tbl %>% 
  filter(median_z > 5 | median_z < -5) %>% 
  filter(percent_samples >= 0.25)
# Bubble plot per type -----------------------------------------------------------
## Custom function to summarize coloc per condition of interest
Make_coloc_per_condition <- function(sample_sets) {
  # Check input
  stopifnot(is.list(sample_sets), !is.null(names(sample_sets)))
  
  # Core function for one condition
  process_one_set <- function(samples, an_condition) {
    pair_totals <- summary_samples_flat %>%
      filter(sample %in% samples) %>%
      distinct(sample, from_cell_type, to_cell_type) %>%
      count(from_cell_type, to_cell_type, name = "n_total")
    
    TME <- c("Mac", "Vasc", "Neuron", "Oligo")
    cancer <- c("cAC", "cMES", "cMESHyp", "cOPC", "cUndiff")
    
    summary_tbl <- summary_samples_flat %>%
      filter(sample %in% samples) %>%
      filter(!from_cell_type %in% c("excluded"),
             !to_cell_type %in% c("excluded")) %>% 
      group_by(from_cell_type, to_cell_type) %>%
      summarise(
        mean_z   = mean(z_score,  na.rm = TRUE),
        median_z = median(z_score, na.rm = TRUE),
        n_sig    = sum(sig,        na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      left_join(pair_totals, by = c("from_cell_type", "to_cell_type")) %>%
      mutate(
        percent_samples = n_sig / n_total,
        pair            = paste0(from_cell_type, "_", to_cell_type),
        condition       = an_condition,
        coloc_type      = case_when(
          (from_cell_type %in% cancer) & (to_cell_type %in% cancer) ~ "cancer_cancer",
          (from_cell_type %in% TME)    & (to_cell_type %in% TME)    ~ "TME_TME",
          TRUE ~ "cancer_TME"
        )
      )
    return(summary_tbl)
  }
  
  # Map over all sample sets and combine
  results <- purrr::imap_dfr(sample_sets, process_one_set)
  
  return(results)
}

## Custom plotting function
Plot_coloc_condition <- function(data = combined_tbl,
                                 facet_order = NULL,
                                 n_rows = 1) {
  
  if (!is.null(facet_order)) {
    data <- data %>%
      mutate(condition = factor(condition, levels = facet_order))
  }
  
  ggplot(data %>% 
           filter(from_cell_type != to_cell_type) %>% 
           drop_na() %>% 
           filter(!from_cell_type %in% c("excluded"),
                  !to_cell_type %in% c("excluded")),
         aes(x = factor(to_cell_type,   levels = order_cell_type4),
             y = factor(from_cell_type, levels = order_cell_type4),
             col = median_z,
             size = percent_samples)) +
    geom_tile(col = "black", size = 0, fill = "white") +
    geom_point(shape = 19) +
    scale_color_gradientn(
      colors = hotmap,
      values = color_values,
      limits = c(-15, 15),
      oob    = scales::squish,
      na.value = "white"
    ) +
    scale_size(range = c(1, 8), limits = c(0.1, 1)) +
    theme_classic() +
    guides(size = guide_legend(override.aes = list(color = "black"))) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, size = "significant in\nfraction of samples") +
    facet_wrap(~condition, nrow = n_rows)
}

## Define and name set of samples
sample_sets <- list(
  IDH_A = samples_IDH_A,
  IDH_O = samples_IDH_O
)

## Run function
combined_tbl <- Make_coloc_per_condition(sample_sets)

## Plot with facets
Plot_coloc_condition()

# Bubble plot grade ---------------------------------------------------------------
## Define and name set of samples
sample_sets <- list(
  low_grade = samples_low,
  mid_grade = samples_mid,
  high_grade = samples_high
)

## Run function
combined_tbl <- Make_coloc_per_condition(sample_sets)

## PLot with facets
Plot_coloc_condition(data = combined_tbl, 
                     facet_order = c("low_grade", "mid_grade", "high_grade"))

# Bubble plot type & grade ----------------------------------------------------------
## Define and name set of samples
sample_sets <- list(
  IDH_A_G2 = samples_IDH_A_G2,
  IDH_A_G3 = samples_IDH_A_G3,
  IDH_A_G4 = samples_IDH_A_G4,
  IDH_O_G2 = samples_IDH_O_G2,
  IDH_O_G3 = samples_IDH_O_G3
)

## Run function
combined_tbl <- Make_coloc_per_condition(sample_sets)

## PLot with facets
Plot_coloc_condition(data = combined_tbl,n_rows = 2)

# Heatmap interaction type per IDH-A grade ------------------------------------------------
an_col_coloc_type <- c(
  "cancer_cancer" = "#ff0066",
  "cancer_TME"    = "#d1aac2",
  "TME_TME"       = "#328c97"
)

## Define and name set of samples
sample_sets <- list(
  IDH_A_G2 = samples_IDH_A_G2,
  IDH_A_G3 = samples_IDH_A_G3,
  IDH_A_G4 = samples_IDH_A_G4
)

# sample_sets <- list(
#   IDH_A = samples_IDH_A,
#   IDH_O = samples_IDH_O
# )


## Run function
combined_tbl <- Make_coloc_per_condition(sample_sets)

filtered_tbl <- combined_tbl %>%
  # 1. Filter by z-score threshold and percent and remove self pairs
  filter(abs(median_z) > 5, percent_samples > 0.5) %>%
  filter(from_cell_type != to_cell_type) %>% 
  mutate(interaction_type = case_when(median_z > 5 ~ "attraction", TRUE ~ "avoidance")) %>% 
  
  # 2. Create a symmetrical key to identify reciprocal pairs
  mutate(pair_key = pmap_chr(list(from_cell_type, to_cell_type), ~paste(sort(c(...)), collapse = "_"))) 

  # 3. For each symmetric pair, keep the row with the most extreme median_z
  # group_by(condition, pair_key) %>%
  # slice_max(order_by = median_z, n = 1, with_ties = FALSE) %>%
  # ungroup()


p1 <- filtered_tbl %>%
  # group_by(condition) %>% distinct(pair_key, .keep_all = TRUE) %>% ungroup() %>% ## If reciproce pairs should be excluded
  ggplot(aes(x = condition, y = pair_key, fill = coloc_type)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = an_col_coloc_type) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_wrap(~interaction_type)

p2 <- filtered_tbl %>% 
  # group_by(condition) %>% distinct(pair_key, .keep_all = TRUE) %>% ungroup() %>% ## If reciproce pairs should be excluded
  count(condition, interaction_type, name = "interaction_count") %>% 
  ggplot(., aes(x = condition, y = interaction_count, group=1)) +
  geom_line() +
  geom_point() +
  facet_wrap(~interaction_type) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) 
  
(p2 / p1) + plot_layout(heights = c(1, 3))


# Heatmap interaction type per all grade ---------------------------------------------------------------------
an_col_coloc_type_reci <- c(
  "cancer_cancer_uni" = "#ff0066",
  "cancer_TME_uni"    = "#d1aac2",
  "TME_TME_uni"       = "#328c97",
  "cancer_cancer_reci" = "#ff0066",
  "cancer_TME_reci"    = "#d1aac2",
  "TME_TME_reci"       = "#328c97"
)

## Note: the interactions in G2 are different between IDH-A and IDH-O. Hence, we calculate the interactions
## seperately per type and then sum them up.

sample_sets <- list(
  IDH_A_G2 = samples_IDH_A_G2,
  IDH_A_G3 = samples_IDH_A_G3,
  IDH_A_G4 = samples_IDH_A_G4,
  IDH_O_G2 = samples_IDH_O_G2,
  IDH_O_G3 = samples_IDH_O_G3
)

## Run function
combined_tbl <- Make_coloc_per_condition(sample_sets)

filtered_tbl <- combined_tbl %>%
  # 1. Filter by z-score threshold and percent and remove self pairs
  filter(abs(median_z) > 5, percent_samples > 0.5) %>%
  filter(from_cell_type != to_cell_type) %>% 
  mutate(interaction_type = case_when(median_z > 5 ~ "attraction", TRUE ~ "avoidance")) %>% 
  mutate(condition = case_when(
    condition %in% c("IDH_A_G2", "IDH_O_G2") ~ "low_grade",
    condition %in% c("IDH_A_G3") ~ "mid_grade",
    condition %in% c("IDH_A_G4", "IDH_O_G3") ~ "high_grade")) %>% 
  mutate(condition = factor(condition, levels = c("low_grade", "mid_grade", "high_grade"))) %>% 
  
  # 2. Create a symmetrical key to identify reciprocal pairs
  mutate(pair_key = pmap_chr(list(from_cell_type, to_cell_type), ~paste(sort(c(...)), collapse = "_"))) %>% 
  group_by(pair_key, condition, coloc_type) %>% ## Add reciprocity
  mutate(
    reciprocity = n_distinct(from_cell_type) > 1  # e.g., Mac→Vasc and Vasc→Mac
  ) %>%
  ungroup() %>% 
  mutate(
    coloc_type_reci = case_when(
      reciprocity ~ paste0(coloc_type, "_reci"),
      TRUE ~ paste0(coloc_type,"_uni")
    ))

p1 <- filtered_tbl %>%
  ggplot(aes(x = condition, y = pair_key, fill = coloc_type_reci)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = an_col_coloc_type_reci) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_wrap(~interaction_type)

p2 <- filtered_tbl %>% 
  count(condition, interaction_type, name = "interaction_count") %>% 
  mutate(condition = factor(condition, levels = c("low_grade", "mid_grade", "high_grade"))) %>% 
  ggplot(., aes(x = condition, y = interaction_count, group=1)) +
  geom_line(linewidht=3) +
  geom_point(size=5) +
  facet_wrap(~interaction_type) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA)
  ) 

(p2 / p1) + plot_layout(heights = c(1, 3))

# Modifications for Figure ------------------------------------------------

##Lineplot count interactions
filtered_tbl %>% 
  count(condition, interaction_type, name = "interaction_count") %>% 
  mutate(
    # rename condition levels
    condition = recode(condition,
                       "low_grade" = "low",
                       "mid_grade" = "mid",
                       "high_grade" = "high"),
    condition = factor(condition, levels = c("low", "mid", "high")),
    # rename facet labels
    interaction_type = recode(interaction_type,
                              "attraction" = "Interaction",
                              "avoidance" = "Repulsion")
  ) %>% 
  ggplot(aes(x = condition, y = interaction_count, group = 1,
             color = interaction_type)) +
  geom_line(linewidth = 2) +
  geom_point(size = 5) +
  facet_wrap(~interaction_type, nrow=2) +
  scale_color_manual(
    values = c("Interaction" = "#b2182b", "Repulsion" = "#2166ac")
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"  # hide redundant legend since facet names explain color
  ) +
  labs(x = NULL, y = "Interaction count")

## Heatmap interactions
ordered_tbl <- filtered_tbl %>%
  group_by(pair_key) %>%
  summarize(order_val = first(coloc_type_reci)) %>%
  arrange(order_val) %>%
  mutate(pair_key = factor(pair_key, levels = unique(pair_key)))

filtered_tbl %>%
  left_join(ordered_tbl %>% select(pair_key, order_val), by = "pair_key") %>%
  mutate(
    pair_key = factor(pair_key, levels = levels(ordered_tbl$pair_key)),
    condition = recode(condition,
                       "low_grade" = "low",
                       "mid_grade" = "mid",
                       "high_grade" = "high"),
    condition = factor(condition, levels = rev(c("low", "mid", "high"))),
    interaction_type = recode(interaction_type,
                              "attraction" = "Interaction",
                              "avoidance" = "Repulsion")
  ) %>%
  ggplot(aes(y = condition, x = pair_key, fill = coloc_type_reci)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = an_col_coloc_type_reci) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_wrap(~interaction_type, nrow=2) +
  labs(y="Cat_grade", x="Interaction pair", fill="Interaction category")

# Final -------------------------------------------------------------------
pal <- c(
  "Cancer-Cancer" = "#ff0066",
  "Cancer-TME"    = "#a5506d",
  "TME-TME"       = "#328c97"
)


collapsed_tbl <- ordered_tbl %>%
  mutate(
    order_val = str_remove(order_val, "_(reci|uni)$"),   # drop suffix
    order_val = case_when(
      order_val == "TME_TME"         ~ "TME-TME",
      order_val == "cancer_cancer"   ~ "Cancer-Cancer",
      order_val == "cancer_TME"      ~ "Cancer-TME",
      TRUE ~ order_val
    )
  )

filtered_tbl %>%
  left_join(collapsed_tbl %>% select(pair_key, order_val), by = "pair_key") %>%
  mutate(
    pair_key = factor(pair_key, levels = levels(ordered_tbl$pair_key)),
    condition = recode(condition,
                       "low_grade" = "Low",
                       "mid_grade" = "Mid",
                       "high_grade" = "High"),
    condition = factor(condition, levels = rev(c("Low", "Mid", "High"))),
    interaction_type = recode(interaction_type,
                              "attraction" = "Interaction",
                              "avoidance" = "Repulsion")
  ) %>%
  ggplot(aes(y = condition, x = pair_key, fill = order_val)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_wrap(~interaction_type, nrow=2) +
  labs(y="Cat. Grade", x="Interaction pair", fill="Interaction category")

## Without repulsions
p1 <- filtered_tbl %>% 
  count(condition, interaction_type, name = "interaction_count") %>% 
  mutate(
    # rename condition levels
    condition = recode(condition,
                       "low_grade" = "Low",
                       "mid_grade" = "Mid",
                       "high_grade" = "High"),
    condition = factor(condition, levels = c("Low", "Mid", "High")),
    # rename facet labels
    interaction_type = recode(interaction_type,
                              "attraction" = "Interaction",
                              "avoidance" = "Repulsion")
  ) %>% 
  filter(interaction_type == "Interaction") %>% 
  ggplot(aes(x = condition, y = interaction_count, group = 1,
             color = interaction_type)) +
  geom_line(linewidth = 2) +
  geom_point(size = 5) +
  scale_color_manual(
    values = c("Interaction" = "#b2182b", "Repulsion" = "#2166ac")
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"  # hide redundant legend since facet names explain color
  ) +
  labs(x = NULL, y = "# Consensus interactions")


p2 <- filtered_tbl %>%
  left_join(collapsed_tbl %>% select(pair_key, order_val), by = "pair_key") %>%
  mutate(
    pair_key = factor(pair_key, levels = levels(ordered_tbl$pair_key)),
    condition = recode(condition,
                       "low_grade" = "Low",
                       "mid_grade" = "Mid",
                       "high_grade" = "High"),
    condition = factor(condition, levels = rev(c("Low", "Mid", "High"))),
    interaction_type = recode(interaction_type,
                              "attraction" = "Interaction",
                              "avoidance" = "Repulsion")
  ) %>%
  filter(interaction_type == "Interaction") %>% 
  ggplot(aes(y = condition, x = pair_key, fill = order_val)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(y="Categortical grade", x="Interaction pair", fill="Interaction category")

(p1+p2) + plot_layout(widths = c(1, 3))

