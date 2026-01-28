# Density & Purity --------------------------------------------------------
## Function
Analyze_density <- function(
    cells, sample_set, title = "Relative Density & Purity (vs mid-zone, binned)",
    thickness = 200,
    malignant_levels = c("cMESHyp","cMES","cAC2","cAC1","cOPC","cUndiff","cNPC"),
    bin_breaks = seq(-2, 2, by = 1),
    bin_labels = paste0("Bin", 1:4),
    reverse_x = TRUE,
    distance_bin = c("junction1","junction2"),
    sample_palette = an_col_lines,
    per_sample_fit = c("gam","line"),
    k_gam = 4,
    per_sample_bin_um = 400,
    drop_no_midzone = TRUE,
    verbose = TRUE,
    ylimits_density = c(-50, 60),   #facet-specific limits for density (%)
    ylimits_purity  = c(-0.3, 0.2)  #facet-specific limits for purity (Δ)
) {
  per_sample_fit <- match.arg(per_sample_fit)
  distance_bin <- match.arg(distance_bin)
  dist_col <- rlang::sym(distance_bin)
  
  bin_centers <- (head(bin_breaks, -1) + tail(bin_breaks, -1))/2
  
  sample_palette <- c("#6BAED6", "#A1D99B", "#F2DF89", "#BDB2FF", "#FC9272", "#FFD8E7", "#DAB894")
  
  ta <- cells %>%
    dplyr::filter(sample %in% sample_set, !is.na(!!dist_col))
  
  roi_areas <- ta %>%
    dplyr::filter(!(cell_type %in% c("excluded"))) %>%
    dplyr::group_by(sample, !!dist_col) %>%
    dplyr::summarise(
      min_x = min(centroid_x, na.rm = TRUE),
      max_x = max(centroid_x, na.rm = TRUE),
      min_y = min(centroid_y, na.rm = TRUE),
      max_y = max(centroid_y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      length_um_x = max_x - min_x,
      length_um_y = max_y - min_y,
      area_mm2 = dplyr::case_when(
        sample == "HBT204A" ~ (length_um_y * thickness) / 1e6,
        TRUE                ~ (length_um_x * thickness) / 1e6
      )
    )
  
  density_per <- ta %>%
    dplyr::filter(!(cell_type %in% c("excluded"))) %>%
    dplyr::group_by(sample, !!dist_col) %>%
    dplyr::summarise(cell_count = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(roi_areas, by = c("sample", distance_bin)) %>%
    dplyr::mutate(cell_density_mm2 = cell_count / area_mm2) %>%
    dplyr::select(sample, !!dist_col, cell_density_mm2)
  
  purity_per <- ta %>%
    dplyr::filter(!(cell_type %in% c("excluded","low"))) %>%
    dplyr::group_by(sample, !!dist_col) %>%
    dplyr::summarise(
      purity = sum(cell_type %in% malignant_levels) / dplyr::n(),
      .groups = "drop"
    )
  
  df_metrics <- dplyr::left_join(density_per, purity_per, by = c("sample", distance_bin)) %>%
    dplyr::mutate(distance_mm = -!!dist_col / 1000)
  
  plot_df <- df_metrics %>%
    tidyr::pivot_longer(
      cols = c("cell_density_mm2","purity"),
      names_to = "metric", values_to = "value"
    )
  
  plot_df_popbins <- plot_df %>%
    dplyr::mutate(
      bin = cut(distance_mm, breaks = bin_breaks, labels = bin_labels,
                include.lowest = TRUE, right = TRUE)
    )
  
  baseline_df <- plot_df_popbins %>%
    dplyr::filter(bin %in% c("Bin2","Bin3")) %>%
    dplyr::group_by(sample, metric) %>%
    dplyr::summarise(baseline = mean(value, na.rm = TRUE), .groups = "drop")
  
  plot_df_scaled <- plot_df_popbins %>%
    dplyr::left_join(baseline_df, by = c("sample","metric")) %>%
    dplyr::mutate(
      value_rel = dplyr::case_when(
        metric == "cell_density_mm2" ~ dplyr::if_else(is.finite(baseline) & baseline > 0,
                                                      100 * (value / baseline - 1),
                                                      NA_real_),
        metric == "purity" ~ value - baseline,
        TRUE ~ NA_real_
      )
    )
  
  if (drop_no_midzone) {
    good_samples <- baseline_df %>%
      dplyr::filter(is.finite(baseline)) %>%
      dplyr::distinct(sample, metric)
    
    missing_info <- dplyr::anti_join(
      plot_df_scaled %>% dplyr::distinct(sample, metric),
      good_samples, by = c("sample","metric")
    )
    
    if (verbose && nrow(missing_info) > 0) {
      message(
        "Dropping ", dplyr::n_distinct(missing_info$sample),
        " sample(s) without mid-zone (Bin2/3) baseline: ",
        paste(unique(missing_info$sample), collapse = ", ")
      )
    }
    
    plot_df_scaled <- plot_df_scaled %>%
      dplyr::semi_join(good_samples, by = c("sample","metric"))
  }
  
  fine_w_mm <- per_sample_bin_um / 1000
  x_range <- range(plot_df_scaled$distance_mm, na.rm = TRUE)
  
  floor_to <- function(x, w) w * floor(x / w)
  ceil_to  <- function(x, w) w * ceiling(x / w)
  x_min <- floor_to(x_range[1], fine_w_mm)
  x_max <- ceil_to(x_range[2],  fine_w_mm)
  
  fine_breaks <- seq(x_min, x_max, by = fine_w_mm)
  if (length(fine_breaks) < 2L) {
    fine_breaks <- c(x_min, x_min + fine_w_mm)
  }
  fine_centers <- (head(fine_breaks, -1) + tail(fine_breaks, -1))/2
  
  plot_df_fine <- plot_df_scaled %>%
    dplyr::mutate(
      fine_bin = cut(distance_mm, breaks = fine_breaks, include.lowest = TRUE, right = TRUE)
    )
  
  sample_fine <- plot_df_fine %>%
    dplyr::group_by(sample, metric, fine_bin) %>%
    dplyr::summarise(
      value_rel = stats::median(value_rel, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      center_x = fine_centers[as.integer(fine_bin)]
    ) %>%
    dplyr::filter(is.finite(center_x), is.finite(value_rel))
  
  bin_summary <- plot_df_scaled %>%
    dplyr::group_by(metric, bin) %>%
    dplyr::summarise(
      mean_value = mean(value_rel, na.rm = TRUE),
      se_value   = stats::sd(value_rel, na.rm = TRUE)/sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bin_num  = as.integer(stringr::str_remove(bin, "Bin")),
      center_x = bin_centers[bin_num]
    )
  
  sample_popbin <- plot_df_scaled %>%
    dplyr::group_by(sample, metric, bin) %>%
    dplyr::summarise(value_rel = stats::median(value_rel, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(bin))
  
  mk_paired <- function(metric_name) {
    sample_popbin %>%
      dplyr::filter(metric == metric_name, bin %in% c("Bin1","Bin4")) %>%
      tidyr::pivot_wider(names_from = bin, values_from = value_rel) %>%
      tidyr::drop_na(Bin1, Bin4)
  }
  
  paired_summary <- function(df) {
    if (nrow(df) < 3) return(NA_real_)
    stats::t.test(df$Bin1, df$Bin4, paired = TRUE)$p.value
  }
  
  majority_consensus <- function(df) {
    if (!all(c("Bin1","Bin4") %in% names(df))) return(NA)
    diffs <- df$Bin4 - df$Bin1
    diffs <- diffs[is.finite(diffs) & diffs != 0]
    if (length(diffs) == 0) return(FALSE)
    prop_same <- max(sum(diffs > 0), sum(diffs < 0)) / length(diffs)
    isTRUE(prop_same > 0.5)
  }
  
  df_cell <- mk_paired("cell_density_mm2")
  df_pur  <- mk_paired("purity")
  
  res_cell <- paired_summary(df_cell)
  res_pur  <- paired_summary(df_pur)
  
  cons_cell <- majority_consensus(df_cell)
  cons_pur  <- majority_consensus(df_pur)
  
  y_max_by_metric <- sample_fine %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(ymax = max(value_rel, na.rm = TRUE), .groups = "drop")
  
  annot_df <- tibble::tibble(
    metric = c("cell_density_mm2","purity"),
    pval   = c(res_cell, res_pur),
    consensus = c(cons_cell, cons_pur)
  ) %>%
    dplyr::mutate(
      label = sprintf(
        "p = %.3f\nConsensus: %s",
        pval,
        ifelse(is.na(consensus), "FALSE", ifelse(consensus, "TRUE", "FALSE"))
      )
    ) %>%
    dplyr::left_join(y_max_by_metric, by = "metric") %>%
    dplyr::mutate(
      x = 0,
      y = ymax * 0.99
    )
  
  # --- Helper: clamp values for plotting to facet-specific limits
  clamp_vec <- function(val, metric) {
    lo <- ifelse(metric == "cell_density_mm2", ylimits_density[1], ylimits_purity[1])
    hi <- ifelse(metric == "cell_density_mm2", ylimits_density[2], ylimits_purity[2])
    pmax(pmin(val, hi), lo)
  }
  
  sample_fine_plot <- sample_fine %>%
    dplyr::mutate(value_plot = clamp_vec(value_rel, metric))
  
  bin_summary_plot <- bin_summary %>%
    dplyr::mutate(
      mean_plot = clamp_vec(mean_value, metric),
      ymin_plot = clamp_vec(mean_value - se_value, metric),
      ymax_plot = clamp_vec(mean_value + se_value, metric)
    )
  
  # Blank anchors to force exact per-facet limits (no padding)
  limits_df <- tibble::tibble(
    metric = rep(c("cell_density_mm2","purity"), each = 2),
    y = c(ylimits_density[1], ylimits_density[2],
          ylimits_purity[1],  ylimits_purity[2]),
    x = 0
  )
  
  p <- ggplot2::ggplot()
  
  if (per_sample_fit == "line") {
    p <- p +
      ggplot2::geom_line(
        data = sample_fine_plot,
        ggplot2::aes(x = center_x, y = value_plot, color = sample, group = interaction(sample, metric)),
        linewidth = 0.6, alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = sample_fine_plot,
        ggplot2::aes(x = center_x, y = value_plot, color = sample),
        size = 1.1, alpha = 0.7
      )
  } else if (per_sample_fit == "gam") {
    p <- p +
      ggplot2::geom_point(
        data = sample_fine_plot,
        ggplot2::aes(x = center_x, y = value_plot, color = sample),
        size = 0.9, alpha = 0.45
      ) +
      ggplot2::geom_smooth(
        data = sample_fine_plot,
        ggplot2::aes(x = center_x, y = value_plot, color = sample, group = interaction(sample, metric)),
        method = "gam",
        formula = y ~ s(x, k = k_gam),
        method.args = list(method = "REML", select = TRUE),
        se = FALSE, linewidth = 0.7, alpha = 0.95
      )
  }
  
  p <- p +
    ggplot2::geom_smooth(
      data = bin_summary_plot,
      ggplot2::aes(x = center_x, y = mean_plot, group = metric),
      method = "loess", span = 0.9, se = FALSE, color = "black", linewidth = 2.5
    ) +
    ggplot2::geom_point(
      data = bin_summary_plot,
      ggplot2::aes(x = center_x, y = mean_plot),
      color = "black", size = 3
    ) +
    ggplot2::geom_errorbar(
      data = bin_summary_plot,
      ggplot2::aes(x = center_x, ymin = ymin_plot, ymax = ymax_plot),
      width = 0.08, color = "black"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_text(
      data = annot_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0.5, vjust = 0, size = 4, lineheight = 0.95, color = "black"
    ) +
    ggplot2::geom_blank(
      data = limits_df,
      ggplot2::aes(x = x, y = y)   # anchors facet-specific limits
    ) +
    ggplot2::facet_wrap(
      ~ metric,
      scales = "free_y",
      labeller = ggplot2::labeller(
        metric = c(
          cell_density_mm2 = "Cell Density",
          purity           = "Purity"
        )
      )
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_color_manual(values = sample_palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.5),
      strip.text   = ggplot2::element_text(size = 14)
    ) +
    ggplot2::labs(
      title = title,
      x = "Distance to Junction (mm)",
      y = NULL,
      color = "Sample"
    )
  
  if (reverse_x) {
    p <- p + ggplot2::scale_x_reverse(limits = c(max(bin_breaks), min(bin_breaks)))
  } else {
    p <- p + ggplot2::scale_x_continuous(limits = c(min(bin_breaks), max(bin_breaks)))
  }
  
  stats_out <- tibble::tibble(
    metric    = c("cell_density_mm2","purity"),
    pval      = c(annot_df$pval[annot_df$metric=="cell_density_mm2"],
                  annot_df$pval[annot_df$metric=="purity"]),
    consensus = c(annot_df$consensus[annot_df$metric=="cell_density_mm2"],
                  annot_df$consensus[annot_df$metric=="purity"])
  )
  
  list(
    plot = p,
    bin_summary = bin_summary,     # summaries on original scaled values
    sample_fine = sample_fine,     # fine-binned per-sample values (unclamped)
    stats = stats_out,
    fine_breaks = fine_breaks,
    ylimits = list(density = ylimits_density, purity = ylimits_purity)
  )
}

tmp1 <- Analyze_density(
  cells,
  samples_WGMJ_IDH_A,
  distance_bin = "junction1",
  reverse_x = F,
  title = NULL,
  per_sample_fit = "line",
  per_sample_bin_um = 1000
)
tmp3 <- Analyze_density(
  cells,
  samples_WWMB_IDH_A,
  distance_bin = "junction2",
  reverse_x = F,
  title = NULL,
  per_sample_fit = "line",
  per_sample_bin_um = 1000
)
tmp1 <- tmp1$plot + labs(x=NULL, y="cells per mm2")
tmp3 <- tmp3$plot + labs(y="cells per mm2") + theme(strip.text = element_blank())

tmp1 / tmp3 +
  plot_layout(
              guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

## Center facets around zero
center_facets_around_zero <- function(p, facet_var = "metric") {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  pp  <- gb$layout$panel_params
  
  get_ylim <- function(ppi) {
    if (!is.null(ppi$y.range)) return(ppi$y.range)
    if (!is.null(ppi$y$range)) return(ppi$y$range)
    if (!is.null(ppi$y$limits)) return(ppi$y$limits)
    c(NA_real_, NA_real_)
  }
  
  lims <- lapply(pp, get_ylim)
  lims_df <- tibble::tibble(
    PANEL = lay$PANEL,
    ymin  = vapply(lims, `[`, numeric(1), 1),
    ymax  = vapply(lims, `[`, numeric(1), 2)
  )
  
  lims_df[[facet_var]] <- lay[[facet_var]]
  lims_df <- lims_df %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::summarise(
      half_range = max(abs(c(ymin, ymax)), na.rm = TRUE),
      .groups = "drop"
    )
  
  limits_long <- lims_df %>%
    tidyr::pivot_longer(cols = "half_range", values_to = "v", names_to = "drop") %>%
    dplyr::select(-drop) %>%
    tidyr::crossing(bound = c("low","high")) %>%
    dplyr::mutate(y = ifelse(bound == "low", -v, v), x = 0)
  
  p + 
    ggplot2::geom_blank(data = limits_long, ggplot2::aes(x = x, y = y)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)))
}

tmp1_centered <- center_facets_around_zero(tmp1, facet_var = "metric")
tmp3_centered <- center_facets_around_zero(tmp3, facet_var = "metric")
(tmp1_centered / tmp3_centered) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

# Cell type composition ---------------------------------------------------
Analyze_cell_comp <- function(
    cells, sample_set,
    title = "Cell-type Δfraction vs mid-zone across junction",
    subtitle = "Black = smoothed Bin means of Δ; dashed red = 0 (mid-zone baseline)",
    thickness = 200,
    levels_cell_type_frac = c("cOPC","Neuron","Oligo","Astro","Vasc","cAC1","cAC2","Mac","cUndiff"),
    exclude_facets = c("Tcell_CD4","Tcell_CD8","neutro","MES2"),
    bin_breaks = seq(-2, 2, by = 1),
    bin_labels = paste0("Bin", 1:4),
    sample_palette = an_col_lines,
    reverse_x = TRUE,
    x_limits = c(-2, 2),
    y_pad = 1.05,
    distance_bin = c("junction1","junction2","CTJI"),
    per_sample_bin_um = 400,
    per_sample_fit = c("gam","line"),
    k_gam = 4,
    drop_no_midzone = TRUE,
    verbose = TRUE
){
  per_sample_fit <- match.arg(per_sample_fit)
  distance_bin <- match.arg(distance_bin); if (identical(distance_bin,"CTJI")) distance_bin <- "junction2"
  dist_col <- rlang::sym(distance_bin)
  
  # Population bin centers
  bin_centers <- (head(bin_breaks,-1) + tail(bin_breaks,-1))/2
  extreme_span_mm <- 3
  
  sample_palette <- c("#6BAED6", "#A1D99B", "#F2DF89", "#BDB2FF", "#FC9272", "#FFD8E7", "#DAB894")
  
  # 1) Cohort filter
  ta <- cells %>%
    dplyr::filter(sample %in% sample_set, !is.na(!!dist_col))
  
  # 2) Per-sample cell-type fractions along distance (absolute fractions)
  cell_type_summary_sample <- ta %>%
    dplyr::group_by(sample, !!dist_col, cell_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups="drop") %>%
    dplyr::group_by(sample, !!dist_col) %>%
    dplyr::mutate(total = sum(count), fraction = count/total) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(cell_type),
                  !(cell_type %in% exclude_facets),
                  cell_type %in% levels_cell_type_frac) %>%
    dplyr::mutate(
      cell_type = factor(cell_type, levels = levels_cell_type_frac),
      distance_mm = -(!!dist_col)/1000
    )
  
  # 3) Assign population bins (Bin1..Bin4) for baseline + summaries
  cell_type_binned <- cell_type_summary_sample %>%
    dplyr::mutate(
      bin = cut(distance_mm, breaks = bin_breaks, labels = bin_labels,
                include.lowest = TRUE, right = TRUE)
    )
  
  # ---- Per-sample mid-zone baseline (Bin2&3) on absolute fractions ----
  baseline_df <- cell_type_binned %>%
    dplyr::filter(bin %in% c("Bin2","Bin3")) %>%
    dplyr::group_by(sample, cell_type) %>%
    dplyr::summarise(baseline = mean(fraction, na.rm = TRUE), .groups = "drop")
  
  # Join baseline and compute Δfraction vs mid-zone (relative space)
  ct_rel <- cell_type_binned %>%
    dplyr::left_join(baseline_df, by = c("sample","cell_type")) %>%
    dplyr::mutate(
      delta_frac = fraction - baseline
    )
  
  # Optionally drop samples without Bin2/3 baseline for a given cell type
  if (drop_no_midzone) {
    keep_pairs <- baseline_df %>% dplyr::distinct(sample, cell_type)
    missing <- dplyr::anti_join(
      ct_rel %>% dplyr::distinct(sample, cell_type), keep_pairs, by = c("sample","cell_type")
    )
    if (verbose && nrow(missing) > 0) {
      message(
        "Dropping ", dplyr::n_distinct(missing$sample),
        " sample(s) lacking mid-zone baseline for at least one cell type: ",
        paste(unique(missing$sample), collapse = ", ")
      )
    }
    ct_rel <- ct_rel %>% dplyr::semi_join(keep_pairs, by = c("sample","cell_type"))
  }
  
  # ---- Fine per-sample binning on a global grid (default 400 µm) ----
  fine_w_mm <- per_sample_bin_um / 1000
  x_range <- range(ct_rel$distance_mm, na.rm = TRUE)
  floor_to <- function(x, w) w * floor(x / w)
  ceil_to  <- function(x, w) w * ceiling(x / w)
  x_min <- floor_to(x_range[1], fine_w_mm)
  x_max <- ceil_to(x_range[2],  fine_w_mm)
  fine_breaks <- seq(x_min, x_max, by = fine_w_mm)
  if (length(fine_breaks) < 2L) fine_breaks <- c(x_min, x_min + fine_w_mm)
  fine_centers <- (head(fine_breaks, -1) + tail(fine_breaks, -1))/2
  
  ct_rel_fine <- ct_rel %>%
    dplyr::mutate(
      fine_bin = cut(distance_mm, breaks = fine_breaks, include.lowest = TRUE, right = TRUE)
    )
  
  # Per-sample-per-cell_type aggregation in fine bins (median is robust)
  sample_fine <- ct_rel_fine %>%
    dplyr::group_by(sample, cell_type, fine_bin) %>%
    dplyr::summarise(
      delta_frac = stats::median(delta_frac, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      center_x = fine_centers[as.integer(fine_bin)]
    ) %>%
    dplyr::filter(is.finite(center_x), is.finite(delta_frac))
  
  # ---- Cohort summaries in population-bin space (Δfraction) ----
  bin_summary <- ct_rel %>%
    dplyr::group_by(cell_type, bin) %>%
    dplyr::summarise(
      mean_value = mean(delta_frac, na.rm = TRUE),
      se_value   = stats::sd(delta_frac, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bin_num  = as.integer(stringr::str_remove(bin,"Bin")),
      center_x = bin_centers[bin_num]
    )
  
  # "Null" line in Δ space is zero
  hline_df <- tibble::tibble(cell_type = levels_cell_type_frac,
                             null_rate = 0) %>%
    dplyr::filter(cell_type %in% unique(ct_rel$cell_type))
  
  # 4) Bin1 vs Bin4 paired across samples — effect per mm (Δ space)
  paired_wide <- ct_rel %>%
    dplyr::filter(bin %in% c("Bin1","Bin4")) %>%
    dplyr::group_by(cell_type, sample, bin) %>%
    dplyr::summarise(val = mean(delta_frac, na.rm = TRUE), .groups="drop") %>%
    tidyr::pivot_wider(names_from = bin, values_from = val) %>%
    tidyr::drop_na(Bin1, Bin4)
  
  paired_summary_tbl <- function(df){
    delta <- df$Bin4 - df$Bin1; n <- length(delta)
    if (n < 3) return(tibble::tibble(effect_per_mm=NA_real_, ci95_lo_per_mm=NA_real_,
                                     ci95_hi_per_mm=NA_real_, p_t_paired=NA_real_,
                                     p_wilcox_paired=NA_real_, n_pairs=n))
    md <- mean(delta); sd_d <- stats::sd(delta); se <- sd_d/sqrt(n); tcrit <- stats::qt(0.975, df=n-1)
    tibble::tibble(effect_per_mm = md/extreme_span_mm,
                   ci95_lo_per_mm = (md - tcrit*se)/extreme_span_mm,
                   ci95_hi_per_mm = (md + tcrit*se)/extreme_span_mm,
                   p_t_paired = stats::t.test(df$Bin1, df$Bin4, paired=TRUE)$p.value,
                   p_wilcox_paired = stats::wilcox.test(df$Bin1, df$Bin4, paired=TRUE)$p.value,
                   n_pairs = n)
  }
  
  stat_table <- paired_wide %>%
    dplyr::group_by(cell_type) %>% dplyr::group_split() %>%
    purrr::map_dfr(~ paired_summary_tbl(.x) %>% dplyr::mutate(cell_type = unique(.x$cell_type))) %>%
    dplyr::arrange(factor(cell_type, levels = levels_cell_type_frac))
  
  # 5) Per-sample binomial vs per-sample null (kept on absolute fractions, as before)
  per_sample_binom <- cell_type_binned %>%
    dplyr::filter(bin %in% c("Bin1","Bin2","Bin3","Bin4")) %>%
    dplyr::group_by(sample, cell_type) %>%
    dplyr::group_modify(~{
      df <- .x
      null_success <- sum(df$count[df$bin %in% c("Bin2","Bin3")], na.rm=TRUE)
      null_total   <- df %>%
        dplyr::filter(bin %in% c("Bin2","Bin3")) %>%
        dplyr::distinct(!!dist_col, .keep_all=TRUE) %>%
        dplyr::summarise(n=sum(total)) %>% dplyr::pull(n)
      null_rate    <- ifelse(null_total>0, null_success/null_total, NA_real_)
      bin1_success <- sum(df$count[df$bin=="Bin1"], na.rm=TRUE)
      bin1_total   <- df %>% dplyr::filter(bin=="Bin1") %>%
        dplyr::distinct(!!dist_col, .keep_all=TRUE) %>%
        dplyr::summarise(n=sum(total)) %>% dplyr::pull(n)
      Bin1_p <- if (!is.na(null_rate) && bin1_total>0) stats::binom.test(bin1_success, bin1_total, p=null_rate)$p.value else NA_real_
      Bin1_dir <- if (is.na(null_rate) || bin1_total==0) NA_real_ else sign(bin1_success/bin1_total - null_rate)
      bin4_success <- sum(df$count[df$bin=="Bin4"], na.rm=TRUE)
      bin4_total   <- df %>% dplyr::filter(bin=="Bin4") %>%
        dplyr::distinct(!!dist_col, .keep_all=TRUE) %>%
        dplyr::summarise(n=sum(total)) %>% dplyr::pull(n)
      Bin4_p <- if (!is.na(null_rate) && bin4_total>0) stats::binom.test(bin4_success, bin4_total, p=null_rate)$p.value else NA_real_
      Bin4_dir <- if (is.na(null_rate) || bin4_total==0) NA_real_ else sign(bin4_success/bin4_total - null_rate)
      tibble::tibble(Bin1_p, Bin1_dir, Bin4_p, Bin4_dir)
    }) %>% dplyr::ungroup()
  
  per_celltype_vote <- per_sample_binom %>%
    dplyr::group_by(cell_type) %>%
    dplyr::reframe({
      vote <- function(dir, p){
        n_tot <- sum(is.finite(p)); n_dir <- sum(!is.na(dir))
        if (n_tot==0 || n_dir==0) return("NS")
        n_up <- sum(dir>0, na.rm=TRUE); n_down <- sum(dir<0, na.rm=TRUE)
        maj  <- dplyr::case_when(n_up > n_dir/2 ~ 1, n_down > n_dir/2 ~ -1, TRUE ~ NA)
        if (is.na(maj)) return("NS")
        frac <- sum(p<0.05 & dir==maj, na.rm=TRUE) / n_tot
        dplyr::case_when(frac >= 0.75 ~ "**", frac > 0.50 ~ "*", TRUE ~ "NS")
      }
      tibble::tibble(Bin1_vote = vote(Bin1_dir, Bin1_p),
                     Bin4_vote = vote(Bin4_dir, Bin4_p))
    })
  
  # 6) Annotations (only votes), positioned in Δ space
  y_annot <- max(sample_fine$delta_frac, na.rm=TRUE); y_annot <- ifelse(is.finite(y_annot), y_annot*y_pad, NA_real_)
  annot_sig_df <- per_celltype_vote %>% dplyr::transmute(cell_type, x=0, y=y_annot, label=paste0(Bin1_vote," ",Bin4_vote))
  
  # 7) Plot (Δ space)
  p <- ggplot2::ggplot()
  
  if (per_sample_fit == "line") {
    p <- p +
      ggplot2::geom_line(
        data = sample_fine,
        ggplot2::aes(x = center_x, y = delta_frac, color = sample, group = interaction(sample, cell_type)),
        linewidth = 0.5, alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = sample_fine,
        ggplot2::aes(x = center_x, y = delta_frac, color = sample),
        size = 1.0, alpha = 0.7
      )
  } else if (per_sample_fit == "gam") {
    p <- p +
      ggplot2::geom_point(
        data = sample_fine,
        ggplot2::aes(x = center_x, y = delta_frac, color = sample),
        size = 0.9, alpha = 0.45
      ) +
      ggplot2::geom_smooth(
        data = sample_fine,
        ggplot2::aes(x = center_x, y = delta_frac, color = sample, group = interaction(sample, cell_type)),
        method = "gam",
        formula = y ~ s(x, k = k_gam),
        method.args = list(method = "REML", select = TRUE),
        se = FALSE, linewidth = 0.7, alpha = 0.95
      )
  }
  
  p <- p +
    ggplot2::geom_smooth(
      data = bin_summary,
      ggplot2::aes(x = center_x, y = mean_value, group = cell_type),
      method = "loess", span = 1, se = FALSE, color = "black", linewidth = 2
    ) +
    ggplot2::geom_point(
      data = bin_summary,
      ggplot2::aes(x = center_x, y = mean_value),
      color = "black", size = 3
    ) +
    ggplot2::geom_errorbar(
      data = bin_summary,
      ggplot2::aes(x = center_x, ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.08, color = "black"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_text(
      data = annot_sig_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 3.2
    ) +
    ggplot2::facet_wrap(~cell_type, nrow = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = sample_palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = ggplot2::element_rect(color="black", fill=NA, size=0.5),
                   strip.text   = element_text(size = 12)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Distance to Junction (mm)",
      y = "Δ Cell type fraction",
      color = "Sample"
    )
  
  if (reverse_x) {
    p <- p + ggplot2::scale_x_reverse(limits = rev(sort(x_limits)))
  } else {
    p <- p + ggplot2::scale_x_continuous(limits = sort(x_limits))
  }
  
  # Return
  list(
    plot = p,
    stat_table = stat_table,          # in Δ space (Bin1 vs Bin4 per cell type)
    per_sample_binom = per_sample_binom,  # as before (absolute-count binomial tests)
    per_celltype_vote = per_celltype_vote,
    bin_summary = bin_summary,        # cohort Δ means by population bin
    hline_df = hline_df,              # 0 baseline for Δ space
    sample_fine = sample_fine,        # fine-binned per-sample Δ values
    fine_breaks = fine_breaks         # the fine grid used
  )
}

tmp1 <- Analyze_cell_comp(
  cells = cells, 
  sample_set = samples_WGMJ_IDH_A,
  title = NULL, 
  levels_cell_type_frac = c("Neuron","Oligo","Astro","Vasc","Mac","cOPC","cAC1","cAC2","cUndiff"),
  per_sample_fit = "line",
  per_sample_bin_um = 800,
  subtitle = NULL,
  reverse_x = F, 
  x_limits = c(-2, 2), 
  y_pad = 0.2)


tmp3 <- Analyze_cell_comp(
  cells = cells, 
  sample_set = samples_WWMB_IDH_A,
  title = NULL,
  levels_cell_type_frac = c("Neuron","Oligo","Astro","Vasc","Mac","cOPC","cAC1","cAC2","cUndiff"),
  distance_bin = "CTJI",
  per_sample_fit = "line",
  per_sample_bin_um = 800,
  subtitle = NULL,
  reverse_x = F, 
  x_limits = c(-2, 2), 
  y_pad = 0.2)

tmp1 <- tmp1$plot +
  labs(x=NULL)
tmp3 <- tmp3$plot
# 
# tmp1 / tmp3 +
#   plot_layout(
#     guides = "collect") &
#   theme(legend.position = "right",
#         legend.box = "vertical")

## Put both plots on the same ylimitis per cell type
# Helper to extract per-facet y-limits from a ggplot with facet_wrap(scales="free_y")
extract_facet_ylims <- function(p, facet_var = "cell_type") {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  pp  <- gb$layout$panel_params
  
  get_ylim <- function(ppi) {
    # be robust across ggplot2 versions
    if (!is.null(ppi$y.range)) return(ppi$y.range)
    if (!is.null(ppi$y$range)) return(ppi$y$range)
    if (!is.null(ppi$y$limits)) return(ppi$y$limits)
    if (!is.null(ppi$y.range_c)) return(ppi$y.range_c)
    c(NA_real_, NA_real_)
  }
  
  lims <- lapply(pp, get_ylim)
  lims_df <- tibble::tibble(
    PANEL = lay$PANEL,
    ymin  = vapply(lims, `[`, numeric(1), 1),
    ymax  = vapply(lims, `[`, numeric(1), 2)
  )
  
  # attach facet label (e.g., cell_type) to each PANEL
  if (!facet_var %in% names(lay)) {
    stop("Facet variable '", facet_var, "' not found in plot layout.")
  }
  lims_df[[facet_var]] <- lay[[facet_var]]
  lims_df %>% select(all_of(facet_var), ymin, ymax) %>% distinct()
}

# 1) Get y-limits per cell_type from tmp1 (AFTER you built tmp1$plot)
ylims_from_tmp1 <- extract_facet_ylims(tmp1, facet_var = "cell_type")

# Optional: small padding so lines don’t touch the edge
ylims_from_tmp1 <- ylims_from_tmp1 %>%
  mutate(pad = pmax(1e-6, 0.02 * (ymax - ymin + 1e-6)),
         ymin = ymin - pad, ymax = ymax + pad) %>%
  select(-pad)

# 2) Turn limits into two invisible points per facet to anchor scales
limits_long <- ylims_from_tmp1 %>%
  pivot_longer(c(ymin, ymax), names_to = "which", values_to = "y") %>%
  mutate(x = 0)   # any x within your x_limits is fine

# 3) Apply those limits to both plots (tmp1 keeps its own look; tmp3 gets aligned)
tmp1_aligned <- tmp1 +
  geom_blank(data = limits_long, aes(x = x, y = y)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = NULL)

tmp3_aligned <- tmp3 +
  geom_blank(data = limits_long, aes(x = x, y = y)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)))

# 4) Patchwork combine
(tmp1_aligned / tmp3_aligned) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

## Center facets around zero
center_facets_around_zero <- function(p, facet_var = "metric") {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  pp  <- gb$layout$panel_params
  
  get_ylim <- function(ppi) {
    if (!is.null(ppi$y.range)) return(ppi$y.range)
    if (!is.null(ppi$y$range)) return(ppi$y$range)
    if (!is.null(ppi$y$limits)) return(ppi$y$limits)
    c(NA_real_, NA_real_)
  }
  
  lims <- lapply(pp, get_ylim)
  lims_df <- tibble::tibble(
    PANEL = lay$PANEL,
    ymin  = vapply(lims, `[`, numeric(1), 1),
    ymax  = vapply(lims, `[`, numeric(1), 2)
  )
  
  lims_df[[facet_var]] <- lay[[facet_var]]
  lims_df <- lims_df %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::summarise(
      half_range = max(abs(c(ymin, ymax)), na.rm = TRUE),
      .groups = "drop"
    )
  
  limits_long <- lims_df %>%
    tidyr::pivot_longer(cols = "half_range", values_to = "v", names_to = "drop") %>%
    dplyr::select(-drop) %>%
    tidyr::crossing(bound = c("low","high")) %>%
    dplyr::mutate(y = ifelse(bound == "low", -v, v), x = 0)
  
  p + 
    ggplot2::geom_blank(data = limits_long, ggplot2::aes(x = x, y = y)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)))
}

tmp1_centered <- center_facets_around_zero(tmp1_aligned, facet_var = "cell_type")
tmp3_centered <- center_facets_around_zero(tmp3_aligned, facet_var = "cell_type")
(tmp1_centered / tmp3_centered) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

# cOPC state changes ------------------------------------------------------
Analyze_marker_means_BinVsNull <- function(
    cells, rmat, sample_set,
    markers = c("SOX9", "OLIG2", "SOX10", "PDGFRA", "SOX2","CD24", "EGFR", "ATRX"),
    cell_type_filter = "cOPC",
    desired_order = markers,
    title = "Marker mean expression across junction (cOPC)",
    subtitle = NULL,
    bin_breaks = seq(-2, 2, by = 1),
    bin_labels = paste0("Bin", 1:4),
    reverse_x = TRUE,
    x_limits = c(-2, 2),
    y_limits = c(-0.4, 0.4),
    annot_top_frac = 0.05,  # fraction of y-range to place text below top
    distance_bin = c("junction1","junction2"),
    an_col_lines = c("#6BAED6", "#A1D99B", "#F2DF89", "#BDB2FF", "#FC9272", "#FFD8E7", "#DAB894"),
    relative = TRUE,
    # NEW: per-sample fine-binning control (um)
    per_sample_bin_um = 400
){
  distance_bin <- match.arg(distance_bin)
  dist_col <- rlang::sym(distance_bin)
  bin_centers <- (head(bin_breaks,-1) + tail(bin_breaks,-1))/2
  extreme_span_mm <- 3
  
  # 1) Filter + 2) expression
  ta <- cells %>% dplyr::filter(sample %in% sample_set, !is.na(!!dist_col))
  expr_mat <- t(rmat[markers, ta %>% dplyr::pull(cell_name), drop = FALSE]); colnames(expr_mat) <- markers
  ta <- cbind(ta, expr_mat) %>% tibble::as_tibble()
  
  # 3) per-sample mean (absolute scale)
  sample_summary <- ta %>%
    dplyr::filter(cell_type == !!cell_type_filter) %>%
    dplyr::select(sample, !!dist_col, dplyr::all_of(markers)) %>%
    tidyr::pivot_longer(dplyr::all_of(markers), names_to = "marker", values_to = "score") %>%
    dplyr::group_by(sample, !!dist_col, marker) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(marker = factor(marker, levels = desired_order),
                  distance_mm = - (!!dist_col) / 1000)
  
  # 4) population bins (Bin1..Bin4)
  marker_binned <- sample_summary %>%
    dplyr::mutate(
      bin = cut(distance_mm, breaks = bin_breaks, labels = bin_labels,
                include.lowest = TRUE, right = TRUE)
    )
  
  bin_summary <- marker_binned %>%
    dplyr::group_by(marker, bin) %>%
    dplyr::summarise(mean_value = mean(mean_score, na.rm = TRUE),
                     se_value   = stats::sd(mean_score, na.rm = TRUE)/sqrt(dplyr::n()),
                     .groups = "drop") %>%
    dplyr::mutate(bin_num = as.integer(stringr::str_remove(bin, "Bin")),
                  center_x = bin_centers[bin_num])
  
  hline_df <- bin_summary %>%
    dplyr::filter(bin %in% c("Bin2","Bin3")) %>%
    dplyr::group_by(marker) %>%
    dplyr::summarise(null_level = mean(mean_value, na.rm = TRUE), .groups = "drop")
  
  # relative scaling vs per-marker mid-zone (Bin2&3)
  if (relative){
    sample_summary_rel <- sample_summary %>% dplyr::left_join(hline_df, by = "marker") %>%
      dplyr::mutate(mean_score_rel = mean_score - null_level)
    bin_summary_rel <- bin_summary %>% dplyr::left_join(hline_df, by = "marker") %>%
      dplyr::mutate(mean_value_rel = mean_value - null_level)
  } else {
    sample_summary_rel <- dplyr::mutate(sample_summary, mean_score_rel = mean_score)
    bin_summary_rel    <- dplyr::mutate(bin_summary,    mean_value_rel = mean_value)
  }
  
  # ---- NEW: Per-sample fine-binning on a global grid (default 400 µm) ----
  fine_w_mm <- per_sample_bin_um / 1000
  x_range <- range(sample_summary_rel$distance_mm, na.rm = TRUE)
  
  floor_to <- function(x, w) w * floor(x / w)
  ceil_to  <- function(x, w) w * ceiling(x / w)
  x_min <- floor_to(x_range[1], fine_w_mm)
  x_max <- ceil_to(x_range[2],  fine_w_mm)
  
  fine_breaks <- seq(x_min, x_max, by = fine_w_mm)
  if (length(fine_breaks) < 2L) fine_breaks <- c(x_min, x_min + fine_w_mm)
  fine_centers <- (head(fine_breaks, -1) + tail(fine_breaks, -1))/2
  
  sample_summary_fine <- sample_summary_rel %>%
    dplyr::mutate(
      fine_bin = cut(distance_mm, breaks = fine_breaks, include.lowest = TRUE, right = TRUE)
    ) %>%
    dplyr::group_by(sample, marker, fine_bin) %>%
    dplyr::summarise(
      mean_score_rel = stats::median(mean_score_rel, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(center_x = fine_centers[as.integer(fine_bin)]) %>%
    dplyr::filter(is.finite(center_x), is.finite(mean_score_rel))
  
  # 5) Bin1 vs Bin4 paired effects (table kept; original absolute scale deltas)
  paired_wide <- marker_binned %>%
    dplyr::filter(bin %in% c("Bin1","Bin4")) %>%
    dplyr::group_by(marker, sample, bin) %>%
    dplyr::summarise(val = mean(mean_score, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = bin, values_from = val) %>%
    tidyr::drop_na(Bin1, Bin4)
  
  paired_summary_tbl <- function(df){
    delta <- df$Bin4 - df$Bin1; n <- length(delta)
    if (n < 3) return(tibble::tibble(effect_per_mm = NA, ci95_lo_per_mm = NA, ci95_hi_per_mm = NA,
                                     p_t_paired = NA, p_wilcox_paired = NA, n_pairs = n))
    md <- mean(delta); sd_d <- stats::sd(delta); se <- sd_d/sqrt(n); tcrit <- stats::qt(0.975, df = n-1)
    tibble::tibble(effect_per_mm = md/extreme_span_mm,
                   ci95_lo_per_mm = (md - tcrit*se)/extreme_span_mm,
                   ci95_hi_per_mm = (md + tcrit*se)/extreme_span_mm,
                   p_t_paired = stats::t.test(df$Bin1, df$Bin4, paired = TRUE)$p.value,
                   p_wilcox_paired = stats::wilcox.test(df$Bin1, df$Bin4, paired = TRUE)$p.value,
                   n_pairs = n)
  }
  
  stat_table <- paired_wide %>%
    dplyr::group_by(marker) %>% dplyr::group_split() %>%
    purrr::map_dfr(~paired_summary_tbl(.x) %>% dplyr::mutate(marker = unique(.x$marker))) %>%
    dplyr::select(marker, effect_per_mm, ci95_lo_per_mm, ci95_hi_per_mm, p_t_paired, p_wilcox_paired, n_pairs) %>%
    dplyr::arrange(factor(marker, levels = desired_order))
  
  # per-sample tests for votes (kept as in original on absolute means vs null)
  per_sample_tests <- marker_binned %>%
    dplyr::filter(bin %in% c("Bin1","Bin2","Bin3","Bin4")) %>%
    dplyr::group_by(sample, marker) %>%
    dplyr::group_modify(~{
      df <- .x
      v_null <- df$mean_score[df$bin %in% c("Bin2","Bin3")]
      v_b1   <- df$mean_score[df$bin == "Bin1"]
      v_b4   <- df$mean_score[df$bin == "Bin4"]
      safe_t <- function(a,b){a<-a[is.finite(a)]; b<-b[is.finite(b)]; if(length(a)>=2&&length(b)>=2) stats::t.test(a,b)$p.value else NA_real_}
      p1 <- safe_t(v_b1, v_null); p4 <- safe_t(v_b4, v_null)
      d1 <- ifelse(length(v_b1)>0&&length(v_null)>0, sign(mean(v_b1,na.rm=TRUE)-mean(v_null,na.rm=TRUE)), NA)
      d4 <- ifelse(length(v_b4)>0&&length(v_null)>0, sign(mean(v_b4,na.rm=TRUE)-mean(v_null,na.rm=TRUE)), NA)
      tibble::tibble(Bin1_p = p1, Bin1_dir = d1, Bin4_p = p4, Bin4_dir = d4)
    }) %>% dplyr::ungroup()
  
  per_marker_vote <- per_sample_tests %>%
    dplyr::group_by(marker) %>%
    dplyr::reframe({
      mk_vote <- function(dir,p){
        n_tot <- sum(is.finite(p)); if(n_tot==0) return("NS")
        n_dir <- sum(!is.na(dir)); if(n_dir==0) return("NS")
        n_up <- sum(dir>0,na.rm=TRUE); n_down <- sum(dir<0,na.rm=TRUE)
        major <- dplyr::case_when(n_up>n_dir/2~1, n_down>n_dir/2~-1, TRUE~NA)
        if (is.na(major)) return("NS")
        n_sig <- sum(p<0.05 & dir==major, na.rm=TRUE); frac <- n_sig/n_tot
        dplyr::case_when(frac>=0.75~"**", frac>0.5~"*", TRUE~"NS")
      }
      tibble::tibble(Bin1_vote=mk_vote(Bin1_dir,Bin1_p),
                     Bin4_vote=mk_vote(Bin4_dir,Bin4_p))
    })
  
  # annotation position based on y_limits
  y_range <- y_limits[2] - y_limits[1]
  y_text  <- y_limits[2] - annot_top_frac * y_range
  annot_sig_df <- per_marker_vote %>%
    dplyr::transmute(marker, x = 0, y = y_text, label = paste0(Bin1_vote," ",Bin4_vote))
  
  # ---- Plot (per-sample lines use fine-binned medians) ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = sample_summary_fine,
      ggplot2::aes(x = center_x, y = mean_score_rel, color = sample, group = interaction(sample, marker)),
      linewidth = 0.4, alpha = 0.6
    ) +
    ggplot2::geom_point(
      data = sample_summary_fine,
      ggplot2::aes(x = center_x, y = mean_score_rel, color = sample),
      size = 1.0, alpha = 0.6
    ) +
    ggplot2::geom_smooth(
      data = bin_summary_rel,
      ggplot2::aes(x = center_x, y = mean_value_rel, group = marker),
      method = "loess", span = 1, se = FALSE, color = "black", linewidth = 2
    ) +
    ggplot2::geom_point(
      data = bin_summary_rel,
      ggplot2::aes(x = center_x, y = mean_value_rel),
      color = "black", size = 3
    ) +
    ggplot2::geom_errorbar(
      data = bin_summary_rel,
      ggplot2::aes(x = center_x, ymin = mean_value_rel - se_value, ymax = mean_value_rel + se_value),
      width = 0.08, color = "black"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_text(
      data = annot_sig_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0.5, vjust = 1, size = 3.2
    ) +
    ggplot2::facet_wrap(~marker, nrow = 1) +
    ggplot2::scale_color_manual(values = an_col_lines) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.5)) +
    ggplot2::labs(
      title = if (relative) paste0(title) else title,
      subtitle = subtitle,
      x = "Distance to Junction (mm)",
      y = if (relative) "Δ mean marker score vs null" else "Mean marker score",
      color = "Sample"
    ) +
    ggplot2::scale_y_continuous(limits = y_limits)
  
  if (reverse_x) {
    p <- p + ggplot2::scale_x_reverse(limits = rev(sort(x_limits)))
  } else {
    p <- p + ggplot2::scale_x_continuous(limits = sort(x_limits))
  }
  
  list(
    plot = p,
    stat_table = stat_table,
    bin_summary = bin_summary,
    per_sample_tests = per_sample_tests,
    per_marker_vote = per_marker_vote,
    hline_df = hline_df,
    # NEW: return fine-binned per-sample data & the fine grid
    sample_fine = sample_summary_fine,
    fine_breaks = fine_breaks
  )
}

tmp1 <- Analyze_marker_means_BinVsNull(
  cells = cells, 
  rmat = rmat, 
  sample_set = samples_WGMJ_IDH_A, 
  markers = c(
    "SOX9",
              "OLIG2",
              "SOX10",
              "SOX2",
              "PDGFRA",
              "CD24",
              "ATRX"
              ),
  title = NULL, 
  subtitle = NULL,
  distance_bin = "junction1",    # or "junction2"
  reverse_x = F,
  x_limits = c(-2, 2),
  y_limits = c(-0.3, 0.3),
  relative = T, 
  per_sample_bin_um = 800
)

tmp3 <- Analyze_marker_means_BinVsNull(
  cells, rmat, samples_WWMB_IDH_A,
  title = NULL,
  markers = c(
    "SOX9",
    "OLIG2",
    "SOX10",
    "SOX2",
    "PDGFRA",
    "CD24",
    "ATRX"
  ),
  distance_bin = "junction2",
  reverse_x = F,
  x_limits = c(-2, 2),
  y_limits = c(-0.3, 0.3),
  relative = T,
  per_sample_bin_um = 800
)

tmp1 <- tmp1$plot +
  labs(x=NULL, y="normalized expression")

tmp3 <- tmp3$plot +
  labs(y="normalized expression")

tmp1 / tmp3 +
  plot_layout(
    guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

# extract_facet_ylims <- function(p, facet_var = "cell_type") {
#   gb <- ggplot_build(p)
#   lay <- gb$layout$layout
#   pp  <- gb$layout$panel_params
#   
#   get_ylim <- function(ppi) {
#     # be robust across ggplot2 versions
#     if (!is.null(ppi$y.range)) return(ppi$y.range)
#     if (!is.null(ppi$y$range)) return(ppi$y$range)
#     if (!is.null(ppi$y$limits)) return(ppi$y$limits)
#     if (!is.null(ppi$y.range_c)) return(ppi$y.range_c)
#     c(NA_real_, NA_real_)
#   }
#   
#   lims <- lapply(pp, get_ylim)
#   lims_df <- tibble::tibble(
#     PANEL = lay$PANEL,
#     ymin  = vapply(lims, `[`, numeric(1), 1),
#     ymax  = vapply(lims, `[`, numeric(1), 2)
#   )
#   
#   # attach facet label (e.g., cell_type) to each PANEL
#   if (!facet_var %in% names(lay)) {
#     stop("Facet variable '", facet_var, "' not found in plot layout.")
#   }
#   lims_df[[facet_var]] <- lay[[facet_var]]
#   lims_df %>% select(all_of(facet_var), ymin, ymax) %>% distinct()
# }
# 
# # 1) Get y-limits per cell_type from tmp1 (AFTER you built tmp1$plot)
# ylims_from_tmp1 <- extract_facet_ylims(tmp1, facet_var = "marker")
# 
# # Optional: small padding so lines don’t touch the edge
# ylims_from_tmp1 <- ylims_from_tmp1 %>%
#   mutate(pad = pmax(1e-6, 0.02 * (ymax - ymin + 1e-6)),
#          ymin = ymin - pad, ymax = ymax + pad) %>%
#   select(-pad)
# 
# # 2) Turn limits into two invisible points per facet to anchor scales
# limits_long <- ylims_from_tmp1 %>%
#   pivot_longer(c(ymin, ymax), names_to = "which", values_to = "y") %>%
#   mutate(x = 0)   # any x within your x_limits is fine
# 
# # 3) Apply those limits to both plots (tmp1 keeps its own look; tmp3 gets aligned)
# tmp1_aligned <- tmp1 +
#   geom_blank(data = limits_long, aes(x = x, y = y)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0))) +
#   labs(x = NULL)
# 
# tmp3_aligned <- tmp3 +
#   geom_blank(data = limits_long, aes(x = x, y = y)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0)))
# 
# # 4) Patchwork combine
# (tmp1_aligned / tmp3_aligned) +
#   patchwork::plot_layout(guides = "collect") &
#   theme(legend.position = "right", legend.box = "vertical")
# 
# ## Center facets around zero
# center_facets_around_zero <- function(p, facet_var = "metric") {
#   gb <- ggplot_build(p)
#   lay <- gb$layout$layout
#   pp  <- gb$layout$panel_params
#   
#   get_ylim <- function(ppi) {
#     if (!is.null(ppi$y.range)) return(ppi$y.range)
#     if (!is.null(ppi$y$range)) return(ppi$y$range)
#     if (!is.null(ppi$y$limits)) return(ppi$y$limits)
#     c(NA_real_, NA_real_)
#   }
#   
#   lims <- lapply(pp, get_ylim)
#   lims_df <- tibble::tibble(
#     PANEL = lay$PANEL,
#     ymin  = vapply(lims, `[`, numeric(1), 1),
#     ymax  = vapply(lims, `[`, numeric(1), 2)
#   )
#   
#   lims_df[[facet_var]] <- lay[[facet_var]]
#   lims_df <- lims_df %>%
#     dplyr::group_by(.data[[facet_var]]) %>%
#     dplyr::summarise(
#       half_range = max(abs(c(ymin, ymax)), na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   limits_long <- lims_df %>%
#     tidyr::pivot_longer(cols = "half_range", values_to = "v", names_to = "drop") %>%
#     dplyr::select(-drop) %>%
#     tidyr::crossing(bound = c("low","high")) %>%
#     dplyr::mutate(y = ifelse(bound == "low", -v, v), x = 0)
#   
#   p + 
#     ggplot2::geom_blank(data = limits_long, ggplot2::aes(x = x, y = y)) +
#     ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)))
# }
# 
# tmp1_centered <- center_facets_around_zero(tmp1_aligned, facet_var = "marker")
# tmp3_centered <- center_facets_around_zero(tmp3_aligned, facet_var = "marker")
# (tmp1_centered / tmp3_centered) +
#   patchwork::plot_layout(guides = "collect") &
#   theme(legend.position = "right", legend.box = "vertical")

# IDH-A inf examples ------------------------------------------------------
tmp1 <- Plot_spatial_map(cells_df = cells,
                         samples = c("HBT204B"), 
                         cell_types = cells$ivygap %>% unique(),
                         cell_type_column = "ivygap",
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "discrete",
                         color_by = "ivygap",
                         color = an_col_ivygap, 
                         ivygap_class = c("ct", "inf", "inf_grey", "MVP", "nec", "none"),
                         legend_title = "IvyGAP",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3, 
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 6, alpha = 1))) 

tmp2 <- Plot_spatial_map(cells_df = cells,
                         samples = c("HBT204B"), 
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "gradient",
                         color_by = "junction1",
                         legend_title = "distance bin\n(mm)",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3,
                         rasterize_points = TRUE,
                         raster_dpi = 300)+
  labs(y=NULL)
  # theme(strip.text = element_blank()
        


tmp3 <- Plot_spatial_map(cells_df = cells,
                         samples = c("HBT204A"), 
                         cell_types = cells$ivygap %>% unique(),
                         cell_type_column = "ivygap",
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "discrete",
                         color_by = "ivygap",
                         color = an_col_ivygap, 
                         ivygap_class = c("ct", "inf", "inf_grey", "MVP", "nec", "none"),
                         legend_title = "IvyGAP",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3, 
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 6, alpha = 1)))

tmp4 <- Plot_spatial_map(cells_df = cells,
                         samples = c("HBT204A"), 
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "gradient",
                         color_by = "junction2",
                         legend_title = "distance bin\n(mm)", 
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3,
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  labs(y=NULL)#+
  # theme(strip.text = element_blank())

(tmp1 + tmp2)/(tmp3+tmp4) +
plot_layout(heights = c(1, 1),   # tmp3 gets 2× vertical space
            guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

# Plot_ROI_Stack ----------------------------------------------------------
Plot_roi_stack <- function(
    cells = cells,
    sample_set = samples_WGMJ_IDH_A,
    reverse_x = TRUE,             # logical or named logical c(S1=TRUE, S2=FALSE)
    thickness = 200,
    color_palette = an_col_simple,
    ggtheme = theme_minimal(),
    levels_cell_type = levels_ct,
    order_legend_cell_type = order_legend,
    roi_col = "junction1"               # <- NEW: name of the junction1-like column to use (e.g., "junction1" or "junction2")
) {
  
  `%ni%` <- Negate(`%in%`)
  
  sample_set <- unique(sample_set)
  
  # reverse_x handling
  if (length(reverse_x) == 1L && is.logical(reverse_x)) {
    reverse_tbl <- tibble::tibble(sample = sample_set, reverse = as.logical(reverse_x))
  } else {
    stopifnot(is.logical(reverse_x), !is.null(names(reverse_x)))
    reverse_tbl <- tibble::tibble(sample = names(reverse_x), reverse = as.logical(reverse_x))
  }
  
  # filter + lock facet order + standardize junction1 column to roi_val
  ta <- cells %>%
    dplyr::filter(sample %in% sample_set) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_set),
      roi_val = .data[[roi_col]]            # <- standardize chosen junction1-like column
    )
  
  malignant_levels <- c("cMESHyp","cMES","cAC2","cAC1","cOPC","cUndiff","cNPC")
  
  # Cell-type level handling (no folding section header)
  # Determine stacking levels (bottom→top). If provided, use user vector; else infer.
  if (!is.null(levels_cell_type)) {
    ct_levels <- levels_cell_type
  } else if (!is.null(names(color_palette)) && length(names(color_palette))) {
    ct_levels <- names(color_palette)
  } else {
    ct_levels <- ta %>% dplyr::distinct(cell_type) %>% dplyr::arrange(cell_type) %>% dplyr::pull(cell_type)
  }
  
  # Legend order: if not provided, use stacking order
  if (is.null(order_legend_cell_type)) {
    legend_levels <- ct_levels
  } else {
    legend_levels <- order_legend_cell_type
  }
  
  # Optional: align palette to levels; silently drop missing, keep existing
  pal_vals <- if (!is.null(names(color_palette))) color_palette[ct_levels] else color_palette
  
  # Fractions per junction1-like bin
  frac_per_ROI <- ta %>%
    dplyr::filter(cell_type %ni% c("low","excluded")) %>%
    dplyr::count(sample, roi_val, cell_type, name = "n") %>%
    dplyr::group_by(sample, roi_val) %>%
    dplyr::mutate(fraction = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fraction_neg = dplyr::if_else(cell_type %in% malignant_levels, fraction, -fraction)) %>%
    tidyr::complete(sample, roi_val, cell_type,
                    fill = list(fraction = 0, fraction_neg = 0, n = 0)) %>%
    dplyr::mutate(
      sample    = factor(sample, levels = sample_set),
      cell_type = factor(cell_type, levels = ct_levels)
    ) %>%
    tidyr::drop_na()
  
  # junction1-like areas (mm^2)
  roi_areas <- ta %>%
    dplyr::filter(cell_type %ni% c("excluded")) %>%
    dplyr::group_by(sample, roi_val) %>%
    dplyr::summarize(
      min_x = min(centroid_x, na.rm = TRUE),
      max_x = max(centroid_x, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      width_um = pmax(0, max_x - min_x),
      area_mm2 = (width_um * thickness) / 1e6,
      sample   = factor(sample, levels = sample_set)
    )
  
  # Density (cells/mm^2)
  density_per_ROI <- ta %>%
    dplyr::filter(cell_type %ni% c("excluded")) %>%
    dplyr::count(sample, roi_val, name = "cell_count") %>%
    dplyr::mutate(sample = factor(sample, levels = sample_set)) %>%
    dplyr::left_join(roi_areas, by = c("sample","roi_val")) %>%
    dplyr::mutate(cell_density_mm2 = cell_count / area_mm2)
  
  # Merge density
  frac_per_ROI <- frac_per_ROI %>%
    dplyr::left_join(density_per_ROI, by = c("sample","roi_val")) %>%
    dplyr::mutate(sample = factor(sample, levels = sample_set))
  
  # Purity per junction1-like bin
  purity_per_ROI <- ta %>%
    dplyr::filter(cell_type %ni% c("low","excluded")) %>%
    dplyr::group_by(sample, roi_val) %>%
    dplyr::summarize(purity = sum(cell_type %in% malignant_levels) / dplyr::n(),
                     .groups = "drop") %>%
    dplyr::mutate(sample = factor(sample, levels = sample_set))
  
  # Direction + units
  frac_per_ROI <- frac_per_ROI %>%
    dplyr::left_join(reverse_tbl, by = "sample") %>%
    dplyr::mutate(
      reverse = dplyr::coalesce(reverse, FALSE),
      ROI_dir = dplyr::if_else(reverse, -roi_val, roi_val),
      ROI_mm  = ROI_dir / 1000,
      sample  = factor(sample, levels = sample_set)
    )
  
  purity_per_ROI <- purity_per_ROI %>%
    dplyr::left_join(reverse_tbl, by = "sample") %>%
    dplyr::mutate(
      reverse = dplyr::coalesce(reverse, FALSE),
      ROI_dir = dplyr::if_else(reverse, -roi_val, roi_val),
      ROI_mm  = ROI_dir / 1000,
      sample  = factor(sample, levels = sample_set)
    )
  
  # Secondary axis scaling
  density_k <- frac_per_ROI$cell_density_mm2 / 1000
  max_frac <- suppressWarnings(max(abs(frac_per_ROI$fraction_neg), na.rm = TRUE))
  max_density_k <- suppressWarnings(max(density_k, na.rm = TRUE))
  scaling_factor <- ifelse(is.finite(max_density_k) && max_density_k > 0,
                           max_frac / max_density_k, 1)
  
  x_breaks <- seq(-2, 2, by = 1)
  y2_breaks <- seq(0, 6, by = 1)
  
  # Plot
  p <- ggplot(frac_per_ROI, aes(x = ROI_mm)) +
    geom_area(
      aes(y = fraction_neg, fill = cell_type),
      stat = "identity", position = "stack", alpha = 1
    ) +
    geom_line(
      aes(y = (cell_density_mm2/1000) * scaling_factor, group = 1, color = "Cell Density"),
      linewidth = 1
    ) +
    geom_line(
      data = purity_per_ROI,
      aes(x = ROI_mm, y = purity, color = "Purity"),
      linewidth = 1
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(
      name = "Metrics",
      values = c("Cell Density" = "black", "Purity" = "magenta")
    ) +
    scale_fill_manual(
      values = pal_vals,
      breaks = legend_levels,
      drop   = TRUE,
      name   = "Cell type"
    ) +
    labs(
      y = "Cell type fraction",
      x = "Distance to junction (mm)"
    ) +
    scale_x_continuous(breaks = x_breaks, limits = c(-2, 2)) +
    scale_y_continuous(
      sec.axis = sec_axis(
        ~ . / scaling_factor,
        name = "Density (cells/mm²)",
        breaks = y2_breaks,
        labels = paste0(y2_breaks, "k")
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.title.y.right = element_text(size = 16, color = "black"),
      axis.text.y.right  = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    ) +
    facet_wrap(~ sample, nrow = 1, scales = "free_x", drop = FALSE) +
    ggtheme +
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 2)
    ) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
  return(p)
}

levels_ct <- c(
  "Mac",
  "TcellCD4",
  "TcellCD8",
  "Bcell",
  "Neutro",
  "Vasc",
  "Astro",
  "Oligo",
  "Neuron",
  "cMESHyp",
  "cMES",
  "cUndiff",
  "cNPC",
  "cAC2",
  "cAC1",
  "cOPC",
  "cGEM"
)

order_legend <-  c(
  "cMESHyp",
  "cMES",
  "cUndiff",
  "cNPC",
  "cAC2",
  "cAC1",
  "cOPC",
  "cGEM",
  "Neuron",
  "Oligo",
  "Astro",
  "Vasc",
  "Neutro",
  "Bcell",
  "TcellCD8",
  "TcellCD4",
  "Mac"
)

## IDH-A - White-Grey Junction
Plot_roi_stack(
  cells = cells,
  sample_set = c(samples_WGMJ_IDH_A), 
  roi_col = "junction1",
  reverse_x = c(
    IDHA19 = T,
    IDHA04 = T,
    IDHA05 = T,
    IDHA06 = T
  ),
  thickness = 200
)

## IDH-A - White-White Junction
Plot_roi_stack(
  cells = cells,
  sample_set = c(samples_WWMB_IDH_A),
  roi_col = "junction2",
  reverse_x = c(
    IDHA20 = T,
    IDHA11 = T,
    IDHA13 = T,
    IDHA04 = T
  ),
  thickness = 200
)
## GBM
Plot_roi_stack(
  cells = gbm_cells,
  sample_set = samples_WGMJ_GBM,
  reverse_x = T,
  thickness = 200
)

# Junction Spatial Maps ---------------------------------------------------
Plot_spatial_map <- function(cell_types = c("Neuron", "Oligo"),
                             cells_df = cells,
                             samples = NULL,             # If NULL, all samples are used
                             facet = TRUE,               # TRUE to facet by sample, FALSE otherwise
                             nrows = 3,
                             color = an_col_simple,      # Named vector of colors for discrete mode
                             cell_type_column = "cell_type",
                             ivygap_class = NULL,        # NULL => don't filter by ivygap unless column exists
                             # Coloring control
                             color_mode = c("discrete", "gradient"),
                             color_by = NULL,            # required for gradient
                             # Gradient parameters
                             gradient_palette = hotmap,
                             gradient_values = seq(0, 1, length.out = length(hotmap)),
                             gradient_limits = c(-2000, 2000),
                             gradient_oob = scales::squish,
                             gradient_labels = scales::label_comma(),
                             legend_title = NULL,
                             # Axes toggle & settings (in mm)
                             fix_axes = TRUE,
                             axis_limits_mm = c(0, 10),
                             axis_break_step = 2.5,
                             # NEW: rasterization controls
                             rasterize_points = FALSE,
                             raster_dpi = 300,
                             point_size = 0.3,
                             ggtheme = theme_minimal()) {
  
  color_mode <- match.arg(color_mode)
  cell_type_sym <- rlang::sym(cell_type_column)
  
  # Base data / filters
  plot_data <- cells_df
  if (color_mode == "discrete") {
    plot_data <- dplyr::filter(plot_data, !!cell_type_sym %in% cell_types)
  }
  if ("ivygap" %in% names(cells_df)) {
    if (is.null(ivygap_class)) ivygap_class <- unique(cells_df$ivygap)
    plot_data <- dplyr::filter(plot_data, .data$ivygap %in% ivygap_class)
  }
  if (!is.null(samples)) {
    plot_data <- dplyr::filter(plot_data, .data$sample %in% samples)
    # ensure facet order follows samples
    plot_data$sample <- factor(plot_data$sample, levels = samples)
  } else {
    # if no samples supplied, preserve original order
    plot_data$sample <- factor(plot_data$sample, levels = unique(plot_data$sample))
  }
  if (nrow(plot_data) == 0) {
    warning("No cells left after filtering. Check cell_types, samples, and ivygap_class.")
    return(ggplot2::ggplot() + ggplot2::ggtitle("No cells to plot"))
  }
  
  # Coordinates: µm → mm
  plot_data <- dplyr::mutate(plot_data,
                             centroid_x_mm = centroid_x / 1000,
                             centroid_y_mm = centroid_y / 1000)
  
  # Prepare the point layer (rasterized or vector)
  make_point_layer <- function() {
    if (isTRUE(rasterize_points)) {
      if (!requireNamespace("ggrastr", quietly = TRUE)) {
        stop("Package 'ggrastr' is required for rasterize_points=TRUE. Install via install.packages('ggrastr').")
      }
      ggrastr::geom_point_rast(size = point_size, alpha = 1, raster.dpi = raster_dpi)
    } else {
      ggplot2::geom_point(size = point_size, alpha = 1)
    }
  }
  
  # Build plot based on mode
  if (color_mode == "discrete") {
    legend_title <- `%||%`(legend_title, "cell type")
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = centroid_x_mm, y = centroid_y_mm, color = as.factor(!!cell_type_sym))
    ) +
      make_point_layer()
    
    # Discrete colors
    if (!is.null(color)) {
      vals <- if (!is.null(names(color))) color[cell_types] else color
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) {
        p <- p + ggplot2::scale_color_manual(
          values = vals,
          breaks = intersect(cell_types, unique(plot_data[[cell_type_column]]))
        )
      } else {
        p <- p + ggplot2::scale_color_discrete()
      }
    } else {
      p <- p + ggplot2::scale_color_discrete()
    }
    
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 6, alpha = 1)))
    
  } else {
    # Gradient mode
    if (is.null(color_by)) stop("In gradient mode, set `color_by` to a numeric column (e.g., 'ROI').")
    if (!color_by %in% names(plot_data)) stop(sprintf("Column '%s' not found.", color_by))
    if (!is.numeric(plot_data[[color_by]])) stop(sprintf("Column '%s' must be numeric.", color_by))
    
    legend_title <- `%||%`(legend_title, color_by)
    
    # Convert color variable µm → mm for legend
    plot_data <- dplyr::mutate(plot_data, .color_mm = .data[[color_by]] / 1000)
    limits_mm <- if (!is.null(gradient_limits)) gradient_limits / 1000 else NULL
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = centroid_x_mm, y = centroid_y_mm, color = .data$.color_mm)
    ) +
      make_point_layer()
    
    # Continuous colorbar
    if (!is.null(gradient_palette)) {
      if (is.null(gradient_values)) {
        gradient_values <- seq(0, 1, length.out = length(gradient_palette))
      }
      p <- p + ggplot2::scale_color_gradientn(
        colours = gradient_palette,
        values  = gradient_values,
        limits  = limits_mm,
        oob     = gradient_oob,
        labels  = gradient_labels,
        name    = legend_title
      )
    } else {
      p <- p + ggplot2::scale_color_gradient(
        limits = limits_mm,
        oob    = gradient_oob,
        labels = gradient_labels,
        name   = legend_title
      )
    }
    
    p <- p + ggplot2::guides(colour = ggplot2::guide_colorbar(barheight = grid::unit(80, "pt")))
  }
  
  # helper: robust breaks function (works with reversed scales)
  .breaks_by <- function(step) {
    force(step)
    function(lims) {
      rng <- sort(range(lims, na.rm = TRUE))
      start <- floor(rng[1] / step) * step
      end   <- ceiling(rng[2] / step) * step
      seq(start, end, by = step)
    }
  }
  
  # Axes: either fixed or auto
  if (isTRUE(fix_axes)) {
    bx <- seq(axis_limits_mm[1], axis_limits_mm[2], by = axis_break_step)
    p <- p +
      ggplot2::scale_x_continuous(breaks = bx, limits = axis_limits_mm) +
      ggplot2::scale_y_reverse(breaks = bx, limits = rev(axis_limits_mm))
  } else {
    p <- p +
      ggplot2::scale_x_continuous(breaks = .breaks_by(axis_break_step)) +
      ggplot2::scale_y_reverse(breaks = .breaks_by(axis_break_step))
  }
  
  # Labels/facets/border (keep themes flexible so you can add theme_classic() outside)
  p <- p +
    ggplot2::labs(x = "Position x (mm)", y = "Position y (mm)", color = if (color_mode == "discrete") legend_title else NULL) +
    ggtheme +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
  if (facet) p <- p + ggplot2::facet_wrap(~sample, nrow = nrows)
  
  return(p)
}

# IDHm --------------------------------------------------------------------
## IDHm
tmp1 <- Plot_spatial_map(cells_df = cells,
                         samples = c(samples_WGMJ_IDH_A), 
                         cell_types = cells$ivygap %>% unique(),
                         cell_type_column = "ivygap",
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "discrete",
                         color_by = "ivygap",
                         color = an_col_ivygap, 
                         ivygap_class = c("ct", "inf", "inf_grey", "MVP", "nec", "none"),
                         legend_title = "IvyGAP",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3, 
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 6, alpha = 1))) +
  labs(x=NULL)

tmp2 <- Plot_spatial_map(cells_df = cells,
                         samples = c(samples_WGMJ_IDH_A), 
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "gradient",
                         color_by = "junction1",
                         legend_title = "distance bin\n(mm)",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3,
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  theme(strip.text = element_blank())

tmp3 <- Plot_roi_stack(
  cells = cells,
  sample_set = c(samples_WGMJ_IDH_A),
  roi_col = "junction1",
    reverse_x = c(
      IDHA19 = T,
      IDHA04 = T,
      IDHA05 = T,
      IDHA06 = T
  ),
  thickness = 200
) +
  guides(
    fill = guide_legend(ncol = 2, byrow = TRUE),
    color = guide_legend(ncol = 2, byrow = TRUE)) #+
# theme(strip.text = element_blank())

(tmp1 / tmp2 / tmp3) +
  plot_layout(heights = c(1, 1, 2),   # tmp3 gets 2× vertical space
              guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

## White-White Matter Border
tmp1 <- Plot_spatial_map(cells_df = cells,
                         samples = c(samples_WWMB_IDH_A), 
                         cell_types = cells$ivygap %>% unique(),
                         cell_type_column = "ivygap",
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "discrete",
                         color_by = "ivygap",
                         color = an_col_ivygap, 
                         ivygap_class = c("ct", "inf", "inf_grey", "MVP", "nec", "none"),
                         legend_title = "IvyGAP",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3, 
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 6, alpha = 1))) +
  labs(x=NULL)

tmp2 <- Plot_spatial_map(cells_df = cells,
                         samples = c(samples_WWMB_IDH_A), 
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "gradient",
                         color_by = "junction2",
                         legend_title = "distance bin\n(mm)",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3,
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  theme(strip.text = element_blank())

tmp3 <- Plot_roi_stack(
  cells = cells,
  sample_set = c(samples_WWMB_IDH_A),
  roi_col = "junction2",
  reverse_x = T,
  thickness = 200
) +
  guides(
    fill = guide_legend(ncol = 2, byrow = TRUE),
    color = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(strip.text = element_blank())

(tmp1 / tmp2 / tmp3) +
  plot_layout(heights = c(1, 1, 2),   # tmp3 gets 2× vertical space
              guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

# GBM ---------------------------------------------------------------------
tmp1 <- Plot_spatial_map(cells_df = gbm_cells,
                         samples = samples_WGMJ_GBM,
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "discrete",
                         color_by = "cell_type", 
                         cell_types = c("Neuron", "Oligo"),
                         legend_title = "Cell Type",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3, 
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 6, alpha = 1))) +
  labs(x=NULL)

tmp2 <- Plot_spatial_map(cells_df = gbm_cells,
                         samples = samples_WGMJ_GBM,
                         facet = TRUE, 
                         nrows = 1,
                         color_mode = "gradient",
                         color_by = "junction1",
                         legend_title = "distance bin\n(mm)",
                         fix_axes = F, axis_limits_mm = c(0,12), axis_break_step = 3,
                         rasterize_points = TRUE,
                         raster_dpi = 300) +
  theme(strip.text = element_blank())

tmp3 <- Plot_roi_stack(
  cells = gbm_cells,
  sample_set = samples_WGMJ_GBM,
  roi_col = "junction1",
  reverse_x = T,
  thickness = 200
) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE),
         color = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(strip.text = element_blank())

(tmp1 / tmp2 / tmp3) +
  plot_layout(heights = c(1, 1, 2),   # tmp3 gets 2× vertical space
              guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical")

