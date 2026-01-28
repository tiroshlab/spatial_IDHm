# Load input --------------------------------------------------------------
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

# Run DEPs ----------------------------------------------------------------
### Set filter if Heatmap should be created for IDH-A vs IDH-O, or both
## Per sample
tmp <- cells %>% 
  # filter(type == "IDH_A") %>% 
  filter(cell_type %ni% c("excluded")) %>% 
  select(cell_name, sample, cell_type)
tmp <- map(tmp$sample %>% unique, ~Run_DEPs_per_sample(., cluster = "cell_type", cell_table = tmp))
tmp <- do.call(cbind, tmp)
tmp <- tmp %>% select(order(colnames(tmp)))
tmp2 <- scale(t(tmp), center = TRUE, scale = TRUE)

## For all samples
tmp3 <- Run_DEPs_hm(matrix = rmat[marker_heatmap,], cluster = "cell_type")
tmp3 <-  scale(t(tmp3), center = TRUE, scale = TRUE)

# Heatmap for "cell_type" ------------------------------------------------
marker_cell_type <- c("CD45", "CD11c", "P2Y12" , "MHCII","CD14", "CD163", "CD206", "CD15",
                   "CD4", "CD3", "CD8", "CD69", "CD19",
                   "CD31", "HSPG2", "aSMA", #"TNC", "CD90", 
                   "NeuN", "ATRX", "SOX10", "OLIG2", "PDGFRA", "EGFR", "SOX2", 
                   "VIM", "GFAP", "S100B" ,"SPARC", "PDPN", "CA9", "Ki67")

## All samples
ComplexHeatmap::Heatmap(
  matrix = tmp3[, c(marker_cell_type)], 
  name = "z-score",
  col = circlize::colorRamp2(breaks = seq(-2, 2, length.out = 11), hotmap),
  row_order = c(
    "Mac",
    "Neutro",
    "TcellCD4",
    "TcellCD8",
    "Bcell", 
    "Vasc", 
    "Neuron",
    "Oligo",
    "cOPC",
    "cUndiff",
    "cAC1",
    "cAC2",
    "Astro",
    "cMES",
    "cMESHyp",
    "excluded",
    "low"
  ),
  column_order = marker_cell_type,
  row_names_side = "left",
  row_dend_side = "right"
)

## Per sample
tmp2_row_order <- sapply(rownames(tmp2), function(x) {
  if (str_detect(x, "^Astro"))       "Astro"      else
    if (str_detect(x, "^Bcell"))       "Bcell"      else
      if (str_detect(x, "^cAC1_"))       "cAC1"       else
        if (str_detect(x, "^cAC2_"))       "cAC2"       else
          if (str_detect(x, "^cMES_"))        "cMES"       else
            if (str_detect(x, "^cMESHyp_"))     "cMESHyp"    else
              if (str_detect(x, "cOPC"))         "cOPC"       else
                if (str_detect(x, "cUndiff"))      "cUndiff"    else
                  if (str_detect(x, "excluded"))     "excluded"   else
                    if (str_detect(x, "low"))          "low"        else
                      if (str_detect(x, "Mac"))          "Mac"        else
                        if (str_detect(x, "Neuron"))       "Neuron"     else
                          if (str_detect(x, "Neutro"))       "Neutro"     else
                            if (str_detect(x, "Oligo"))        "Oligo"      else
                              if (str_detect(x, "^TcellCD4_"))   "TcellCD4"  else
                                if (str_detect(x, "^TcellCD8_"))   "TcellCD8"   else
                                  if (str_detect(x, "Vasc"))         "Vasc"       else
                                    NA_character_
})

## Check distribution
table(tmp2_row_order)

Heatmap(tmp2[, marker_cell_type],
        col = circlize::colorRamp2(breaks = seq(-1.5, 1.5, length.out = 11), hotmap),
        row_names_side = "left",
        show_row_names = F,
        row_dend_side = NULL,
        column_order = marker_cell_type,
        row_gap = unit(1, "mm"),
        border = T,
        name = "z-score",
        cluster_row_slices = F,
        cluster_column_slices = F,
        row_split = factor(tmp2_row_order, 
                           levels = c("Mac",
                                      "Neutro",
                                      "TcellCD4",
                                      "TcellCD8",
                                      "Bcell", 
                                      "Vasc", 
                                      "Neuron",
                                      "Oligo",
                                      "cOPC",
                                      "cUndiff",
                                      "cAC1",
                                      "cAC2",
                                      "Astro",
                                      "cMES",
                                      "cMESHyp",
                                      "low")),  # levels sorted alphabetically (or as you desire)
        left_annotation = rowAnnotation(cell_type = tmp2_row_order,
                                        col = list(cell_type = an_col_simple),
                                        show_legend = F),
        column_split = factor(c(
          rep("myeloid", 8),
          rep("lymph", 5),
          rep("vasc", 3),
          rep("brain", 13),
          rep("func", 1)
        ), levels = c("myeloid", "lymph", "vasc", "brain", "func")),
        row_title_rot = 0,
        column_title_gp = gpar(fontsize = 13, fontface = "bold", col = "black"),
        row_title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
        column_names_gp = gpar(col = "black", fontsize = 13, fontfamily = "Helvetica"),
        show_row_dend = F,
        heatmap_legend_param = list(title_gp = gpar(col = "black"),
                                    labels_gp = gpar(col = "black"),
                                    border = "black")
)

# Heatmap for "cell_type2" ------------------------------------------------
marker_cell_type2 <- c("CD45", "CD11c", "P2Y12" , "MHCII","CD14", "CD163", "CD206", "CD15",
                   "CD4", "CD3", "CD8", "CD69", "CD19",
                   "CD31", "HSPG2", "aSMA", "TNC", "CD90", 
                   "NeuN", "ATRX", "SOX10", "OLIG2", "SOX2", "EGFR", 
                   "VIM", "GFAP", "S100B" ,"SPARC", "PDPN", "CA9", "Ki67")

tmp <- cells %>% filter(cell_type2 %ni% "excluded") %>% filter(type == "IDH_O")
tmp <- map(tmp$sample %>% unique, ~Run_DEPs_per_sample(sample = ., cluster = "cell_type2", cell_table = tmp))
tmp <- do.call(cbind, tmp)
tmp <- tmp %>% select(order(colnames(tmp)))
tmp2 <- scale(t(tmp), center = TRUE, scale = TRUE)

## For all samples
tmp3 <- Run_DEPs_hm(matrix = rmat[marker_heatmap,], cluster = "cell_type2")
tmp3 <-  scale(t(tmp3), center = TRUE, scale = TRUE)

## All samples
ComplexHeatmap::Heatmap(
  matrix = tmp3[, c(marker_cell_type2)], 
  name = "z-score",
  col = circlize::colorRamp2(breaks = seq(-2, 2, length.out = 11), hotmap),
  row_order = c(
    "InfMg",
    "GAM",
    "MacScav",
    "MacBorder",
    "Neutro",
    "TcellCD4",
    "TcellCD8",
    "Bcell", 
    "VascBBB",
    "Pericyte",
    "VascAng",
    "Neuron",
    "Oligo",
    "cOPC",
    "cUndiff",
    "cAC1",
    "cAC2",
    "Astro",
    "cMES",
    "cMESHyp",
    "excluded",
    "low"
  ),
  column_order = marker_cell_type2,
  row_names_side = "left",
  row_dend_side = "right"
)

## Per sample
tmp2_row_order <- sapply(rownames(tmp2), function(x) {
  if (str_detect(x, "^InfMg_"))         "InfMg" else
    if (str_detect(x, "^GAM"))           "GAM" else
      if (str_detect(x, "^MacScav_"))      "MacScav" else
        if (str_detect(x, "^MacBorder_"))     "MacBorder" else
          if (str_detect(x, "^Neutro_"))        "Neutro" else
            if (str_detect(x, "^TcellCD4_"))     "TcellCD4" else
              if (str_detect(x, "^TcellCD8_"))     "TcellCD8" else
                if (str_detect(x, "^Bcell_"))         "Bcell" else
                  if (str_detect(x, "^VascBBB_"))          "VascBBB" else
                    if (str_detect(x, "^Pericyte_"))      "Pericyte" else
                      if (str_detect(x, "^VascAng_"))       "VascAng" else
                        if (str_detect(x, "^Neuron_"))        "Neuron" else
                          if (str_detect(x, "^Oligo_"))         "Oligo" else
                            if (str_detect(x, "^cOPC_"))           "cOPC" else
                              if (str_detect(x, "cUndiff"))            "cUndiff" else
                                if (str_detect(x, "^cAC1_"))          "cAC1" else
                                  if (str_detect(x, "^cAC2_"))          "cAC2" else
                                    if (str_detect(x, "^Astro_"))         "Astro" else
                                      if (str_detect(x, "^cMES_"))          "cMES" else
                                        if (str_detect(x, "^cMESHyp_"))          "cMESHyp" else
                                          if (str_detect(x, "low"))             "low" else
                                            NA_character_
})

Heatmap(tmp2[, marker_cell_type2],
        col = circlize::colorRamp2(breaks = seq(-2, 2, length.out = 11), hotmap),
        row_names_side = "left",
        show_row_names = F,
        row_dend_side = NULL,
        column_order = marker_cell_type2,
        row_gap = unit(1, "mm"),
        border = T,
        name = "z-score",
        cluster_row_slices = F,
        cluster_column_slices = F,
        row_split = factor(tmp2_row_order,
                           levels = c("InfMg",
                                      "GAM",
                                      "MacScav",
                                      "MacBorder",
                                      "Neutro",
                                      "TcellCD4",
                                      "TcellCD8",
                                      "Bcell", 
                                      "VascBBB",
                                      "Pericyte",
                                      "VascAng",
                                      "Neuron",
                                      "Oligo",
                                      "cOPC",
                                      "cUndiff",
                                      "cAC1",
                                      "cAC2",
                                      "Astro",
                                      "cMES",
                                      "cMESHyp",
                                      "low")),
        left_annotation = rowAnnotation(cell_type_2 = tmp2_row_order,
                                        col = list(cell_type_2 = an_col_all),
                                        show_legend = F),
        column_split = factor(c(rep("myeloid", 8),
                                rep("lymph", 5),
                                rep("vasc", 5),
                                rep("brain", 12),
                                rep("func", 1)),
                              levels = c("myeloid", "lymph", "vasc", "brain", "func", "ECM", "other")),
        row_title_rot = 0,
        column_title_gp = gpar(fontsize = 13, fontface = "bold", col = "black"),
        row_title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
        column_names_gp = gpar(col = "black", fontsize = 13, fontfamily = "Helvetica"),
        show_row_dend = F,
        heatmap_legend_param = list(title_gp = gpar(col = "black"),
                                    labels_gp = gpar(col = "black"),
                                    border = "black"))
