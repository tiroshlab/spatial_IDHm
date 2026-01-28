# Load libraries, custom functions and colors -------

setwd("/home/projects/tirosh/rouvenho/IDHm/IDHm_github/")
#setwd(".../IDHm_github/")

stopifnot(file.exists("inputs"))

library(ComplexHeatmap)
library(MoMAColors)
library(kohonen)
library(FNN)
library(Rcpp)
library(ggdark)
library(ggrastr)
library(scales)
library(treemapify)
library(ggridges)
# library(ggh4x)
library(ggrastr)
library(tidyverse)
library(arrow)
library(rhdf5)
library(patchwork)

`%ni%` <- Negate(`%in%`)

Open <- function(file = "2_Colors_Markers.R") {
  file.edit({{file}})
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

Count_cell_type <- function(cell_type = "cell_type") {
  cells %>%
    count(!!sym(cell_type)) %>%
    arrange(desc(n)) %>%
    print(n = Inf)
}

hotmap <- c(
  "#053061",
  "#2166ac",
  "#4393c3",
  "#92c5de",
  "#d1e5f0",
  "#f7f7f7",
  "#fddbc7",
  "#f4a582",
  "#d6604d",
  "#b2182b",
  "#67001f"
)

Run_DEPs_hm <- function(matrix = rmat, cell_table = cells, cluster = "pg_cluster") {
  
  # Initialize DEPs_hm as a data frame to preserve row names
  DEPs_hm <- data.frame(row.names = rownames(matrix))
  
  # Iterate over each unique cluster
  for (i in sort(unique(cell_table %>% pull(!!sym(cluster))))) {
    
    # Calculate the difference between means for the current cluster vs. all others
    top_DEP_value <- rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) == i)]) - 
      rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) != i)])
    
    # Add the DEP values as a new column to DEPs_hm
    DEPs_hm[, i] <- top_DEP_value
  }
  
  return(DEPs_hm)
}

Run_DEPs <- function(matrix = rmat, cell_table = cells, cluster = "pg_cluster") {
  # Initialize DEPs as a data frame with the correct number of rows and an empty first column
  DEPs <- data.frame(matrix(nrow = nrow(matrix), ncol = 0, dimnames = list(rownames(matrix), NULL)))
  
  # Iterate over each unique cluster
  for (i in sort(unique(cell_table %>% pull(!!sym(cluster))))) {    
    # Calculate the difference between means for the current cluster vs. all others
    top_DEP_values <- rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) == i)]) - 
      rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) != i)])
    
    # Sort the DEPs by value (decreasing) and get the indices
    top_DEP_ix <- order(top_DEP_values, decreasing = TRUE)
    
    # Select the top DEPs and their corresponding gene names
    DEPs <- cbind(DEPs, Gene = rownames(matrix)[top_DEP_ix], Value = top_DEP_values[top_DEP_ix])
  }
  
  # Assign column names to the DEPs data frame
  cluster_names <- sort(unique(cell_table %>% pull(!!sym(cluster))))
  colnames(DEPs) <- unlist(lapply(cluster_names, function(i) c(paste(i, "Gene"), paste(i, "Value"))))
  
  # Return the DEPs data frame
  return(DEPs)
}

Run_means <- function(matrix = rmat, cell_table = cells, cluster = "pg_cluster") {
  # Initialize DEPs as a data frame with the correct number of rows and row names
  DEPs <- data.frame(matrix(nrow = nrow(matrix), ncol = 0, dimnames = list(rownames(matrix), NULL)))
  
  # Iterate over each unique cluster
  for (i in sort(unique(cell_table %>% pull(!!sym(cluster))))) {    
    # Calculate the row means for the current cluster
    top_DEP_values <- rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) == i)])
    
    # Sort the means by value (decreasing) and get the indices
    top_DEP_ix <- order(top_DEP_values, decreasing = TRUE)
    
    # Select the top means and their corresponding gene names
    DEPs <- cbind(DEPs, Gene = rownames(matrix)[top_DEP_ix], Value = top_DEP_values[top_DEP_ix])
  }
  
  # Assign column names to the DEPs data frame
  cluster_names <- sort(unique(cell_table %>% pull(!!sym(cluster))))
  colnames(DEPs) <- unlist(lapply(cluster_names, function(i) c(paste(i, "Gene"), paste(i, "Value"))))
  
  # Return the DEPs data frame
  return(DEPs)
}

Run_means_hm <- function(matrix = rmat, cell_table = cells, cluster = "pg_cluster") {
  
  # Initialize DEPs_hm as a data frame to preserve row names
  DEPs_hm <- data.frame(row.names = rownames(matrix))
  
  # Iterate over each unique cluster
  for (i in sort(unique(cell_table %>% pull(!!sym(cluster))))) {
    
    # Calculate the difference between means for the current cluster vs. all others
    top_DEP_value <- rowMeans(matrix[, which(cell_table %>% pull(!!sym(cluster)) == i)])
    
    # Add the DEP values as a new column to DEPs_hm
    DEPs_hm[, i] <- top_DEP_value
  }
  
  return(DEPs_hm)
}

Plot_hierarchical_clust <- function(matrix = rmat,
                                    marker = marker_mac,
                                    cluster = "gating",
                                    cell_types = c("Mac_CD206", "Inf_Mac", "Mic", "Mac"),
                                    slice_size = 2000,
                                    color_min = 0,
                                    color_max = 0.9,
                                    n_clusters = 6,
                                    cell_table = cells) {
  
  cell_table <- cells
  
  # Step 1: Subset and subsample the matrix
  tmp <- matrix[marker, cell_table %>% filter(!!sym(cluster) %in% cell_types) %>% pull(cell_name)]
  tmp <- tmp[, sample(colnames(tmp), size = slice_size, replace = FALSE)] # Subsample
  
  # Step 2: Create the column annotation
  an_sample <- HeatmapAnnotation(an_sample = cell_table %>% filter(cell_name %in% colnames(tmp)) %>% pull(sample),
                                 which = "column",
                                 annotation_label = "sample",
                                 annotation_legend_param = list(border = "black"),
                                 border = TRUE,
                                 col = list(an_sample = an_cols_sample)
  )
  
  # Step 3: Generate the main heatmap object
  ht <- Heatmap(tmp,
                name = "score",
                col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                row_names_side = "left",
                row_dend_side = "left",
                show_column_names = FALSE,
                clustering_method_rows = "ward.D",
                clustering_method_columns = "ward.D",
                column_title_side = "bottom",
                column_split = n_clusters, 
                border = TRUE,
                top_annotation = an_sample)
  
  # Step 4: Get the column indices for each group defined by column_split
  column_groups <- column_order(ht)
  
  # Step 5: Compute the average marker expression for each group
  avg_expression <- sapply(column_groups, function(cols) {
    rowMeans(tmp[, cols, drop = FALSE])
  })
  
  colnames(avg_expression) <- as.character(1:dim(avg_expression)[2])
  
  # Step 6: Create the heatmap for the averaged data
  avg_heatmap <- Heatmap(avg_expression,
                         name = "avg expression",
                         col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                         cluster_columns = FALSE,  # Do not cluster columns since they represent groups
                         cluster_rows = FALSE,     # Optional: keep rows unclustered
                         show_column_names = TRUE,
                         column_names_rot = 0,    # Rotate column names by 90 degrees
                         show_row_names = TRUE,
                         width = unit(4, "cm"),
                         border = TRUE)
  
  # Step 7: Combine the two heatmaps
  ht + avg_heatmap
}

Plot_marker_scatter <- function(cell_table = cells,
                                mat = rmat, 
                                cluster = gating,
                                cell_type = "Mac",
                                marker_1 = "CD14", 
                                marker_2 = "CD11c") {
  # Subset the matrix and cells
  cells_filtered <- cell_table %>% 
    filter({{cluster}} %in% c({{cell_type}})) %>% 
    pull(cell_name)
  
  mat_subset <- mat[c(marker_1, marker_2), cells_filtered] %>% 
    t() %>% 
    as.data.frame()
  
  # Add the sample column to the subset data frame
  mat_subset$sample <- cell_table %>% 
    filter({{cluster}} %in% c({{cell_type}})) %>% 
    pull(sample)
  
  # Plotting using tidy evaluation
  ggplot(mat_subset, aes(x = .data[[marker_1]], y = .data[[marker_2]], color = sample)) +
    geom_point(size = 0.5, alpha = 0.5) +
    #geom_hline(yintercept = 0.5, color = "red") +
    labs(color = "Sample") +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    dark_theme_gray()
}

Add_to_cells <- function(category = "immune", cluster = "cell_type") {
  
  subclust_data <- subclust[[category]][[1]]  # Extract the data frame from the list
  
  updated_cells <- cells %>% 
    left_join(subclust_data %>%
                select(cell_name, cluster_column = cell_type),
              by = "cell_name") %>%
    mutate(cluster_column = if_else(is.na(cluster_column), 
                                    true = !!sym(cluster), 
                                    false = cluster_column)) %>%
    select(-!!sym(cluster)) %>%
    rename(!!sym(cluster) := cluster_column)
  
  # Assign the updated tibble back to the global environment
  assign("cells", updated_cells, envir = .GlobalEnv)
}

Export_qpan <- function(cluster = "pg_cluster") {
  qpan <- cells %>% 
    select(image_id = sample,
           cell_id,
           pg_cluster = !!sym(cluster)) # Dynamically select the cluster column
  
  # Create the filename using the first sample name
  file_name <- paste0("qpan_", cells$sample[1], ".csv")
  
  # Write the CSV file
  write_csv(x = qpan, file = file_name, col_names = TRUE)
}

Make_hclust <- function(method = "average", 
                        matrix = DEPs_hm) 
{
  M_dist <- dist(matrix %>% t())
  M_hc  <- hclust(M_dist, method = method)             
  M_hc  <- reorder(as.dendrogram(M_hc), colMeans(matrix %>% t()))
}

Run_DEPs_per_sample <- function(sample, cluster = "gating", matrix = rmat, cell_table = cells) {
  matrix_sample <- matrix[, cell_table %>% 
                            filter(sample %in% c({{sample}})) %>% 
                            filter(!is.na({{cluster}})) %>%
                            pull(cell_name),
                          drop = FALSE]
  
  cells_sample <- cell_table %>% 
    filter(sample %in% c({{sample}})) %>%  
    filter(!is.na({{cluster}}))
  
  DEPs_hm_sample <- tibble(dim(matrix_sample)[1]:dim(matrix_sample)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cells_sample[[cluster]]))) {
    # i iterates all clusters
    in_cluster <- which(i == cells_sample[[cluster]])
    out_cluster <- which(i != cells_sample[[cluster]])
    
    top_DEP_value <- rowMeans(matrix_sample[, in_cluster, drop = FALSE]) - 
      rowMeans(matrix_sample[, out_cluster, drop = FALSE])
    DEPs_hm_sample <- cbind(DEPs_hm_sample, top_DEP_value) # Select top gene names
  }
  
  DEPs_hm_sample <- DEPs_hm_sample[, -1] # Delete redundant first column
  colnames(DEPs_hm_sample) <- paste0(sort(unique(cells_sample[[cluster]])), "_", sample)
  
  return(DEPs_hm_sample)
}

Plot_spatial_map <- function(cell_types = c("Neuron", "Oligo"),
                             cells = cells,
                             samples = NULL,          # If NULL, all samples are used
                             facet = TRUE,            # TRUE to facet by sample, FALSE otherwise
                             nrows = 3,
                             color = an_col_simple,   # Named vector of colors for the cell types (optional)
                             cell_type_column = "cell_type",
                             ivygap_class = NULL) {   # NULL => don't filter by ivygap unless column exists
  
  # Symbol for tidy-eval
  cell_type_sym <- rlang::sym(cell_type_column)
  
  # Base filter on requested cell types
  plot_data <- cells %>%
    dplyr::filter(!!cell_type_sym %in% cell_types)
  
  # Conditionally filter by ivygap if that column exists
  if ("ivygap" %in% names(cells)) {
    if (is.null(ivygap_class)) {
      ivygap_class <- unique(cells$ivygap)
    }
    plot_data <- plot_data %>% dplyr::filter(.data$ivygap %in% ivygap_class)
  }
  
  # Optional sample filter
  if (!is.null(samples)) {
    plot_data <- plot_data %>% dplyr::filter(.data$sample %in% samples)
  }
  
  # Early exit if nothing to plot
  if (nrow(plot_data) == 0) {
    warning("No cells left after filtering. Check cell_types, samples, and ivygap_class.")
    return(ggplot() + ggtitle("No cells to plot"))
  }
  
  # Build plot
  p <- ggplot(plot_data, aes(x = centroid_x, y = centroid_y, color = as.factor(!!cell_type_sym))) +
    geom_point(size = 0.3, alpha = 1) +
    scale_y_reverse() +
    labs(color = "cell type")
  
  # Color handling: use provided colors if possible, else fall back
  if (!is.null(color)) {
    # If named, align to requested cell_types; otherwise use as-is
    vals <- if (!is.null(names(color))) color[cell_types] else color
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      p <- p + scale_color_manual(values = vals,
                                  breaks = intersect(cell_types, unique(plot_data[[cell_type_column]])))
    } else {
      p <- p + scale_color_discrete()
    }
  } else {
    p <- p + scale_color_discrete()
  }
  
  # Legend appearance & theme (keep your original intent)
  p <- p +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    dark_theme_gray()
  
  # Optional faceting
  if (facet) {
    p <- p + facet_wrap(~sample, nrow = nrows)
  }
  
  return(p)
}


Plot_subclust_hm <- function(matrix = rmat,
                             cluster = "gating",
                             markers = marker_mac,
                             cell_type_state = "Mac",
                             color_min = -2,
                             color_max = 2) {
  
  # Filter cells based on the specified cluster and cell_type_state
  selected_cells <- cells %>% 
    dplyr::filter(!!sym(cluster) %in% c(cell_type_state)) %>% 
    dplyr::pull(cell_name)
  
  print(paste("subsetting", length(selected_cells), "cells of type", cell_type_state, sep = " "))
  
  tmp_mat <- matrix[markers, selected_cells] %>% 
    as.matrix() %>% 
    t()
  
  print(paste("dimension of matrix is", dim(tmp_mat)[1], "x", dim(tmp_mat)[2], sep = " "))
  
  #Run_means_hm for sub-clusters
  tmp_means_hm <- Run_means_hm(matrix = tmp_mat %>% t(), cell_table = subclust$Mac$df, cluster = "hc_clust")
  colnames(tmp_means_hm) <- c(1:dim(tmp_means_hm)[2])
  
  tmp_hm_means <- ComplexHeatmap::Heatmap(tmp_means_hm %>% as.matrix(), 
                                          name = "score",
                                          col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                                          column_km = 0,
                                          row_km = 0,
                                          row_names_side = "left",
                                          row_dend_side = "left",
                                          column_names_rot = 0)
  
  # Return the results as a list
  return(list(df = tmp_means_hm,
              hm_means = tmp_hm_means
  )
  )
}

Print_last_n_columns <- function(df = cells, n = 5) {
  df %>% select((ncol(.) - n + 1):ncol(.))
}

# Load inputs -------------------------------------------------------------
## Cell dfs
cells      <- read_parquet("inputs/cells_df.parquet")
gcells     <- read_parquet("inputs/gcells_df.parquet")
gbm_cells  <- read_parquet("inputs/gbm_cells_df.parquet")

## Expression matrices
### Mean intensities
h5ls("inputs/CODEX_IDHm_mean_intensity.h5")

amat <- h5read(
  "inputs/CODEX_IDHm_mean_intensity.h5",
  "CODEX_IDHm_mean_intensity"
)

rn <- h5read(
  "inputs/CODEX_IDHm_mean_intensity.h5",
  "rownames"
)

cn <- h5read(
  "inputs/CODEX_IDHm_mean_intensity.h5",
  "colnames"
)

rownames(amat) <- as.character(rn)
colnames(amat) <- as.character(cn)

## Nimbus Scores IDHm
h5ls("inputs/CODEX_IDHm_nimbus_scores.h5")

rmat <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores.h5",
  "CODEX_IDHm_nimbus_scores"
)

rn <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores.h5",
  "rownames"
)

cn <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores.h5",
  "colnames"
)

rownames(rmat) <- as.character(rn)
colnames(rmat) <- as.character(cn)

## Nimbus Scores GBM
h5ls("inputs/CODEX_IDHm_nimbus_scores_gbm.h5")

gbm_rmat <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores_gbm.h5",
  "CODEX_IDHm_nimbus_scores_gbm"
)

rn <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores_gbm.h5",
  "rownames"
)

cn <- h5read(
  "inputs/CODEX_IDHm_nimbus_scores_gbm.h5",
  "colnames"
)

rownames(gbm_rmat) <- as.character(rn)
colnames(gbm_rmat) <- as.character(cn)

## TmnApp cells
TmnApp_CODEX_cells <- read_parquet("inputs/TmnApp_CODEX_cells.parquet")

# Colors ------------------------------------------------------------------
an_col_simple <- c(
  "cGEM" = "#795c32",
  "Astro" = "#D4D915",
  "Oligo" = "#6e948c",
  "Neuron" = "#251714",
  "cAC1" = "#c43e3b",
  "cAC2" = "#ef7673",
  "cOPC" = "#0066CC",
  "cUndiff" = "#b0799a",
  "cNPC" = "#b0799a",
  "cMESHyp" = "#FFB224",
  "cMES" = "#FFF05A",
  "Vasc" = "#CF1C90",
  "Mac" = "#5e9432",
  "Neutro" = "white",
  "Bcell" = "green",
  "TcellCD8" = "red",
  "TcellCD4" = "red",
  "low"       = "grey"
)
an_col_all <-  c(
  "cUndiff"   = "#78adb7",
  "Neuron"    = "#661f66",
  "cAC1"      = "#c43e3b",
  "cAC2"      = "#ef7673",
  "cGEM" = "#795c32",
  "Astro"     = "#D4D915",
  "cOPC"      = "#0066CC",
  "Oligo"     = "#C993A2",
  "cMESHyp"   = "#ee9b43",
  "cMES"      = "#f3d567",
  "Mac"       = "#4c9a77",
  # endothelial
  "VascBBB"   = "#CF1C90",
  "Vasc"   = "#CF1C90",
  "VascAng"   = "#DA70D6",  
  "Pericyte"  = "#63e6e6",
  # immune
  "GAM"       = "#80BA5A",
  "MacScav"   = "#ff4d6f",
  "MacBorder" = "#f9c000",
  "InfMg"     = "#579ea4",
  "Neutro"    = "white",
  "Bcell"     = "green",
  "TcellCD8"  = "#E41A1C",  
  "TcellCD4"  = "#E41A1C",  
  "low"       = "grey"
)
an_col_mac <- c(
  "GAM" = "#86ad34",
  "MacScav" = "#ff4d6f",
  "MacBorder" = "#df7713",
  "InfMg" = "#579ea4"
  # "TcellCD8" = "red",
  # "TcellCD4" = "red"
)

an_col_vasc <- c(
  "VascBBB" = "#CF1C90",
  "VascAng" = "red",
  "Pericyte" = "#63e6e6"
)

an_col_MP <-  c(
  "cUndiff" = "#78adb7",
  "Neuron" = "#661f66",
  "cAC" = "#e74b47",
  "cOPC" = "#0066CC",
  "Oligo" = "#C993A2",
  "cOXPHOS" = "#033366",
  "cMESAst" = "#ee9b43",
  "cMESHyp" = "#f3d567",
  "InfMac" = "#4c9a77",
  "VascAng" = "#762b35",
  "VascIMEC" = "#af4646",
  "TAM" = "#1b7975",
  "NA" = "white"
)

an_col_MP2 <-
  c(
    "Undiff" = "#78adb7",
    "Neuron" = "#661f66",
    "AC" = "#e74b47",
    "OPC" = "#005187",
    "oligo" = "#C993A2",
    "OXPHOS" = "#172767",
    "MES1" = "#ee9b43",
    "MES2" = "#f3d567",
    "InfMac" = "#4c9a77",
    "VascAng" = "#762b35",
    "VascIMEC" = "#af4646",
    "TAM" = "#1b7975",
    "NA" = "white"
  )

an_col_Tcell <- c(
  "TcellCD4rest" = "#A4CA9A",
  "TcellCD4act" = "#538E9C",
  "TcellCD4ex" = "#213361",
  "TcellCD8rest" = "#FF9898",
  "TcellCD8act" = "#A91E45",
  "TcellCD8ex" = "#431424"
)

an_col_ivygap <- c(
  "artifact" = "black",
  "ct" = "#2A9D3D",
  "inf" = "#d1aac2",
  "inf_grey" = "#a5506d",
  "MVP" = "#ff0066",
  "nec" = "#328c97",
  "none" = "grey"
)

an_cols_sample <- c(
  "IDHO01" = "#1f77b4",      # blue
  "IDHA01" = "#aec7e8",      # light blue
  "IDHA02" = "#ff7f0e",      # orange
  "IDHA03" = "#ffbb78",      # light orange
  "IDHO02" = "#2ca02c",      # green
  "IDHO09" = "#98df8a",    # light green
  "IDHA04" = "#ff9896",    # light red
  "IDHA05" = "#9467bd",    # purple
  "IDHA20" = "#c5b0d5",     # light purple
  "IDHA06" = "#8c564b",     # brown
  "IDHA07" = "#c49c94",     # light brown
  "IDHA08" = "#e377c2",    # pink
  "IDHA09" = "#f7b6d2",    # light pink
  "IDHA18" = "#bcbd22",    # olive
  "IDHA19" = "#c7c7c7",    # light grey
  "IDHA10" = "#7f7f7f",    # grey
  "IDHA11" = "#dbdb8d",    # light olive
  "IDHO10" = "#9edae5",    # light teal
  "IDHO03" = "#393b79",    # dark blue
  "IDHA12" = "#5254a3",    # medium blue
  "IDHA13" = "#6b6ecf",    # light blue-purple
  "IDHO04" = "#9c9ede",     # pale blue-purple
  "IDHO08" = "#8c6d31",     # dark brown
  "IDHA14" = "#bd9e39",      # golden brown
  "IDHO07" = "#e7ba52",      # yellow-brown
  "IDHO05" = "#637939",     # dark olive
  "IDHO06" = "#8ca252",     # olive green
  "IDHA24" = "#b5cf6b",   # light olive green
  "IDHA23" = "#cedb9c",   # pale green
  "IDHA21" = "#ad494a"    # dark red
)

# Marker sets -------------------------------------------------------------
marker_hq <-
  c(
    ## Immune
    "CD11c",
    "CD14",
    "CD15",
    "CD163",
    #"CD19",
    "CD206",
    "CD3",
    "CD45",
    "CD4",
    #"CD69",
    "CD8",
    "CD83",
    #"CD279",
    #"CXCR4",
    "MHCII",
    "P2Y12",
    ## Vasc
    "aSMA",
    "CD31",
    "HSPG2",
    ## Normal
    "GFAP",
    "NeuN",
    "VIM",
    ## Cancer
    "ATRX",
    "EGFR",
    "SOX2",
    "SOX9",
    ## OPC
    "APOD",
    "OLIG2",
    "PDGFRA",
    "SOX10",
    ## MES
    "CA9",
    "GLUT1",
    "NDRG1",
    "PDPN",
    "CHI3L1",
    ## NPC
    "CD24",
    "DCX",
    #"DCXv2",
    "SOX4",
    #"STMN1",
    "TUBB3",
    ## AC
    "APOE",
    "AQP4",
    "PTPRZ1",
    "S100B",
    "SPARC",
    ## Functional
    "IL1B",
    "Ki67",
    ## ECM
    #"ACAN",
    "BCAN",
    "CD44",
    "CD90",
    "COL1",
    "CSPG5",
    "FLT1",
    "FN1",
    "GAP43",
    "HPLN1",
    "MAP2",
    "MT1MMP",
    #"NCAN",
    #"NG2",
    "PLP1",
    "POSTN",
    "SNAP25",
    "TNC",
    "TNR"
    #"VCAN"
  )

marker_gating <-
  c(
    ## Immune
    "CD11c",
    "CD14",
    "CD15",
    "CD163",
    "CD19",
    "CD206",
    "CD3",
    "CD45",
    "CD4",
    "CD69",
    "CD8",
    #"CD83",
    #"CD279",
    #"CXCR4",
    "MHCII",
    "P2Y12",
    ## Vasc
    "aSMA",
    "CD31",
    "HSPG2",
    ## Normal
    #"GFAP",
    "NeuN",
    #"VIM",
    ## Cancer
    "ATRX",
    "EGFR",
    "SOX2",
    #"SOX9",
    ## OPC
    #"APOD",
    "OLIG2",
    #"PDGFRA",
    "SOX10",
    ## MES
    #"CA9",
    #"GLUT1",
    #"NDRG1",
    #"PDPN",
    #"CHI3L1",
    ## NPC
    #"CD24",
    #"DCX",
    #"DCXv2",
    #"SOX4",
    #"STMN1",
    #"TUBB3",
    ## AC
    #"APOE",
    #"AQP4",
    #"PTPRZ1",
    "S100B",
    "SPARC"
    ## Functional
    #"IL1B",
    #"Ki67",
    ## ECM
    # "ACAN",
    # "BCAN",
    # "CD44",
    # "CD90",
    # "COL1",
    # #"CSPG5",
    # "FLT1",
    # "FN1",
    # "GAP43",
    # "HPLN1",
    # "MAP2",
    # "MT1MMP",
    # "NCAN",
    # #"NG2",
    # "PLP1",
    # "POSTN",
    # "SNAP25",
    # "TNC",
    # "TNR",
    # "VCAN"
  )

marker_nimbus <-
  c(
    ## Immune
    "CD11c",
    "CD14",
    "CD15",
    "CD163",
    #"CD19",
    "CD206",
    "CD3",
    "CD45",
    "CD4",
    #"CD69",
    "CD8",
    #"CD83",
    #"CD279",
    #"CXCR4",
    "MHCII",
    "P2Y12",
    ## Vasc
    "aSMA",
    "CD31",
    "HSPG2",
    ## Normal
    "GFAP",
    "NeuN",
    "VIM",
    ## Cancer
    "ATRX",
    "EGFR",
    "SOX2",
    "SOX9",
    ## OPC
    "APOD",
    "OLIG2",
    "PDGFRA",
    "SOX10",
    ## MES
    "CA9",
    "GLUT1",
    "NDRG1",
    "PDPN",
    "CHI3L1",
    ## NPC
    "CD24",
    "DCX",
    #"DCXv2",
    "SOX4",
    #"STMN1",
    "TUBB3",
    ## AC
    "APOE",
    "AQP4",
    "PTPRZ1",
    "S100B",
    "SPARC",
    ## Functional
    "IL1B",
    "Ki67",
    ## ECM
    #"ACAN",
    "BCAN",
    "CD44",
    "CD90",
    "COL1",
    "CSPG5",
    "FLT1",
    "FN1",
    "GAP43",
    "HPLN1",
    "MAP2",
    "MT1MMP",
    #"NCAN",
    #"NG2",
    "PLP1",
    "POSTN",
    "SNAP25",
    "TNC",
    "TNR"
    #"VCAN"
  )

marker_heatmap <-
  c(
    ## Immune
    "CD45",
    "CD11c",
    "MHCII",
    "CD14",
    "CD163",
    "P2Y12",
    "CD206",
    "CD15",
    "CD83",
    "CD19",
    "CD4",
    "CD3",
    "CD69",
    "CD8",
    #"CXCR4",
    ## Vasc
    "CD31",
    "HSPG2",
    "FN1",
    "aSMA",
    "COL1",
    "CD90",
    ## brain
    "NeuN",
    "MAP2",
    "SNAP25",
    "ATRX",
    "SOX10",
    "APOD",
    "OLIG2",
    "SOX2",
    "SOX9",
    "EGFR",
    "SPARC",
    "S100B",
    "GFAP",
    "APOE",
    "AQP4",
    "CD24",
    "DCX",
    #"DCXv2",
    "SOX4",
    #"STMN1",
    "TUBB3",
    "VIM",
    "CHI3L1",
    "GAP43",
    "PDPN",
    "GLUT1",
    "NDRG1",
    "CA9",
    ## Functional
    "Ki67",
    "IL1B",
    ## ECM
    "BCAN",
    "MT1MMP",
    "HPLN1",
    "POSTN",
    "ACAN",
    "VCAN",
    "CSPG5",
    "FLT1",
    #"NCAN",
    #"NG2",
    "TNC",
    "TNR",
    ## others
    "CD44",
    "PLP1",
    "PDGFRA",
    "PTPRZ1",
    "CD279"
  )

marker_heatmap_clean <-
  c(
    ## Immune
    "CD45",
    "CD11c",
    "MHCII",
    "P2Y12",
    "CD14",
    "CD163",
    "CD206",
    "CD15",
    "CD83",
    "CD19",
    "CD4",
    "CD3",
    "CD69",
    "CD8",
    #"CXCR4",
    ## Vasc
    "CD31",
    "HSPG2",
    "GLUT1",
    "FN1",
    "aSMA",
    "COL1",
    "CD90",
    ## brain
    "NeuN",
    "MAP2",
    "SNAP25",
    "ATRX",
    "SOX10",
    "APOD",
    "OLIG2",
    "SOX2",
    "SOX9",
    "EGFR",
    "S100B",
    "SPARC",
    "GFAP",
    "VIM",
    "AQP4",
    "CD44",
    #"DCX",
    #"DCXv2",
    #"SOX4",
    #"STMN1",
    "TUBB3",
    "PDPN",
    "CA9",
    "NDRG1",
    ## Functional
    "Ki67",
    "IL1B",
    ## ECM
    "BCAN",
    "MT1MMP",
    "HPLN1",
    "POSTN",
    "ACAN",
    "VCAN",
    "CSPG5",
    "FLT1",
    #"NCAN",
    #"NG2",
    "TNC",
    "TNR",
    ## others
    "CD24",
    "CHI3L1",
    "GAP43",
    "PLP1",
    "PDGFRA",
    "APOE",
    "PTPRZ1",
    "CD279"
  )

marker_all <-
  c(
    ## Immune
    "CD11c",
    "CD14",
    "CD15",
    "CD163",
    "CD19",
    "CD206",
    "CD3",
    "CD45",
    "CD4",
    "CD69",
    "CD8",
    "CD83",
    "CD279",
    "CXCR4",
    "MHCII",
    "P2Y12",
    ## Vasc
    "aSMA",
    "CD31",
    "HSPG2",
    ## Normal
    "GFAP",
    "NeuN",
    "VIM",
    ## Cancer
    "ATRX",
    "EGFR",
    "SOX2",
    "SOX9",
    ## OPC
    "APOD",
    "OLIG2",
    "PDGFRA",
    "SOX10",
    ## MES
    "CA9",
    "GLUT1",
    "NDRG1",
    "PDPN",
    "CHI3L1",
    ## NPC
    "CD24",
    "DCX",
    #"DCXv2",
    "SOX4",
    #"STMN1",
    "TUBB3",
    ## AC
    "APOE",
    "AQP4",
    "PTPRZ1",
    "S100B",
    "SPARC",
    ## Functional
    "IL1B",
    "Ki67",
    ## ECM
    "ACAN",
    "BCAN",
    "CD44",
    "CD90",
    "COL1",
    "CSPG5",
    "FLT1",
    "FN1",
    "GAP43",
    "HPLN1",
    "MAP2",
    "MT1MMP",
    "NCAN",
    "NG2",
    "PLP1",
    "POSTN",
    "SNAP25",
    "TNC",
    "TNR",
    "VCAN"
  )

marker_simple <- c(
  "CD45",
  "CD11c",
  "P2Y12" ,
  "MHCII",
  "CD14",
  "CD163",
  "CD206",
  "CD15",
  "CD4",
  "CD3",
  "CD8",
  "CD69",
  "CD19",
  "CD31",
  "HSPG2",
  "aSMA",
  "TNC",
  "CD90",
  "NeuN",
  "ATRX",
  "SOX10",
  "OLIG2",
  "SOX2",
  "EGFR",
  "VIM",
  "GFAP",
  "S100B",
  "SPARC",
  "PDPN",
  "CA9",
  "Ki67"
)

# Levels ------------------------------------------------------------------
sample_metadata <- readxl::read_excel(
  path  = "inputs/TableS1_SampleMetadata_v13.xlsx",
  sheet = "This_cohort"
)

samples_IDH_A <- sample_metadata %>% 
  filter(GliomaType %in% c("IDH_A")) %>% 
  pull(SampleID) %>% 
  unique()

samples_IDH_O <- sample_metadata %>% 
  filter(GliomaType %in% c("IDH_O")) %>% 
  pull(SampleID) %>% 
  unique()

samples_WGMJ_IDH_A <- c("IDHA19",
                        "IDHA05",
                        "IDHA06",
                        "IDHA04"
                        )

samples_GGMB_IDH_O <- c("IDHO08", "IDHO05", "IDHO10")

samples_GWMJ_IDH_O <- c("IDHO02", "IDHO06")

samples_WWMB_IDH_A <- c("IDHA04", "IDHA20", "IDHA11", "IDHA13")

samples_WGMJ_GBM <- c(
  "GBM27",
  "GBM28",
  "GBM29",
  "GBM30",
  "GBM31",
  "GBM32",
  "GBM21"
)

samples_IDH_A_G2 <- cells %>% filter(type == "IDH_A", section_grade == 2) %>% pull(sample) %>% unique()
samples_IDH_A_G3 <- cells %>% filter(type == "IDH_A", section_grade == 3) %>% pull(sample) %>% unique()
samples_IDH_A_G4 <- cells %>% filter(type == "IDH_A", section_grade == 4) %>% pull(sample) %>% unique()
samples_IDH_O_G2 <- cells %>% filter(type == "IDH_O", section_grade == 2) %>% pull(sample) %>% unique()
samples_IDH_O_G3 <- cells %>% filter(type == "IDH_O", section_grade == 3) %>% pull(sample) %>% unique()

samples_low <- cells %>% filter(section_grade == 2) %>% pull(sample) %>% unique()
samples_mid <- cells %>% filter(type == "IDH_A", section_grade == 3) %>% pull(sample) %>% unique()
samples_high <- cells %>% 
  filter((type == "IDH_A" & section_grade == 4) | (type == "IDH_O" & section_grade == 3)) %>% 
  pull(sample) %>% unique()
