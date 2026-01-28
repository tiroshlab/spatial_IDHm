# IDH-A -------------------------------------------------------------------
# Myeloid -----------------------------------------------------------------
## Myeloid
keep_cts <- c("InfMg", "GAM", "MacScav", "MacBorder")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type2, type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type2 %in% keep_cts) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type2 = factor(cell_type2, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"mid",
    section_grade == 4 ~"high"
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))

p1 <- ggplot(frac_focus %>% filter(type=="IDH_A"),
       aes(x = cat_grade,
           y = fraction,
           fill = cell_type2)) +
  geom_col(position = "stack") +
  facet_grid(~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_mac,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    x= "Cat_grade",
    y    = "Relative cell-type abundance",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

## Significance IDH_A only
# 1) Per-sample fractions for your chosen cell types & IDH_A only
frac_sample_keep <- cells %>%
  filter(cell_type2 %in% keep_cts,
         type == "IDH_A") %>%
  count(sample, section_grade, cell_type2) %>%       # count cells per sample×section_grade×cell_type2
  group_by(sample) %>%                       
  mutate(frac = n / sum(n)) %>%               # fraction within each sample
  ungroup()

# 2) Pairwise Wilcoxon (2 vs 3, 3 vs 4, AND 2 vs 4) per cell_type2
pairwise_stats_keep <- frac_sample_keep %>%
  group_by(cell_type2) %>%
  do({
    d <- .
    # 2 vs 3
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(2,3)))
    ) %>% mutate(comparison = "2_vs_3")
    # 3 vs 4
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(3,4)))
    ) %>% mutate(comparison = "3_vs_4")
    # 2 vs 4
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(2,4)))
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")  # BH‐correction across all tests
  ) %>%
  arrange(cell_type2, comparison)

# 3) Inspect your results
pairwise_stats_keep %>% arrange(p.value) %>% print(n=Inf)


# T cell ------------------------------------------------------------------
keep_cts <- c("TcellCD4rest", "TcellCD4act", "TcellCD4ex",
              "TcellCD8rest", "TcellCD8act", "TcellCD8ex")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type_Tcell, type, cell_type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type_Tcell %ni% c("other")) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type_Tcell = factor(cell_type_Tcell, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"mid",
    section_grade == 4 ~"high"
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))


p2 <- ggplot(frac_focus %>% filter(type=="IDH_A"),
       aes(x = cat_grade,
           y = fraction,
           fill = cell_type_Tcell)) +
  geom_col(position = "stack") +
  facet_grid(~ cell_type + type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_Tcell,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    y    = "Relative cell-type abundance",
    x    = "Cat_grade",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

##Significance CD4
keep_cd4 <- c("TcellCD4rest", "TcellCD4act", "TcellCD4ex")

# 1) per-sample fractions for CD4
frac_cd4 <- cells %>%
  filter(cell_type_Tcell %in% keep_cd4) %>%
  filter(type == "IDH_A") %>% 
  count(sample, section_grade, type, cell_type_Tcell) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# 2) pairwise Wilcoxon per CD4 subset, approximate p-values
stats_cd4 <- frac_cd4 %>%
  group_by(cell_type_Tcell) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(cell_type_Tcell, comparison)

stats_cd4 %>% arrange(p.value)

## Significance CD8
keep_cd8 <- c("TcellCD8rest", "TcellCD8act", "TcellCD8ex")

# 1) per-sample fractions for CD8 subsets
frac_cd8 <- cells %>%
  filter(cell_type_Tcell %in% keep_cd8) %>%
  filter(type == "IDH_A") %>% 
  count(sample, section_grade, type, cell_type_Tcell) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# 2) pairwise Wilcoxon per CD8 subset, approximate p-values
stats_cd8 <- frac_cd8 %>%
  group_by(cell_type_Tcell) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(cell_type_Tcell, comparison)

# View CD8 results
stats_cd8 %>% arrange(p.value)


# Vasc --------------------------------------------------------------------
keep_cts <- c("VascBBB", "Pericyte", "VascAng")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type2, type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type2 %in% keep_cts) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type2 = factor(cell_type2, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"mid",
    section_grade == 4 ~"high"
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))

p3 <- ggplot(frac_focus %>% filter(type=="IDH_A"),
       aes(x = cat_grade,
           y = fraction,
           fill = cell_type2)) +
  geom_col(position = "stack") +
  facet_grid(~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_vasc,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    y    = "Relative cell-type abundance",
    x    = "Cat_grade",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

## Significance
keep_cts <- c("VascBBB", "Pericyte", "VascAng")

# 1) Per-sample fractions for your vascular subsets (IDH_A only)
frac_vasc <- cells %>%
  filter(cell_type2 %in% keep_cts,
         type == "IDH_A") %>%
  count(sample, section_grade, cell_type2) %>%      # count per sample×section_grade×cell_type2
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%              # fraction within each sample
  ungroup()

# 2) Pairwise Wilcoxon (2 vs 3, 3 vs 4, 2 vs 4) per vascular subset
stats_vasc <- frac_vasc %>%
  group_by(cell_type2) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%  # BH-correction across all tests
  arrange(cell_type2, comparison)

# 3) View results
stats_vasc %>% arrange(p.value)
# Plot all ----------------------------------------------------------------
p1 <- p1 + labs(title = "Myeloid")
p2 <- p2 + labs(title = "T cell")
p3 <- p3 + labs(title = "Vasc")

(p1 + p2 + p3) + plot_layout(widths = c(1, 2, 1)) &
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# IDH-O -------------------------------------------------------------------
# Myeloid -----------------------------------------------------------------
## Myeloid
keep_cts <- c("InfMg", "GAM", "MacScav", "MacBorder")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type2, type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type2 %in% keep_cts) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type2 = factor(cell_type2, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"high",
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))

p1 <- ggplot(frac_focus %>% filter(type=="IDH_O"),
             aes(x = cat_grade,
                 y = fraction,
                 fill = cell_type2)) +
  geom_col(position = "stack") +
  facet_grid(~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_mac,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    x= "Cat_grade",
    y    = "Relative cell-type abundance",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

## Significance IDH_A only
# 1) Per-sample fractions for your chosen cell types & IDH_A only
frac_sample_keep <- cells %>%
  filter(cell_type2 %in% keep_cts,
         type == "IDH_A") %>%
  count(sample, section_grade, cell_type2) %>%       # count cells per sample×section_grade×cell_type2
  group_by(sample) %>%                       
  mutate(frac = n / sum(n)) %>%               # fraction within each sample
  ungroup()

# 2) Pairwise Wilcoxon (2 vs 3, 3 vs 4, AND 2 vs 4) per cell_type2
pairwise_stats_keep <- frac_sample_keep %>%
  group_by(cell_type2) %>%
  do({
    d <- .
    # 2 vs 3
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(2,3)))
    ) %>% mutate(comparison = "2_vs_3")
    # 3 vs 4
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(3,4)))
    ) %>% mutate(comparison = "3_vs_4")
    # 2 vs 4
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade, data = filter(d, section_grade %in% c(2,4)))
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")  # BH‐correction across all tests
  ) %>%
  arrange(cell_type2, comparison)

# 3) Inspect your results
pairwise_stats_keep %>% arrange(p.value) %>% print(n=Inf)


# T cell ------------------------------------------------------------------
keep_cts <- c("TcellCD4rest", "TcellCD4act", "TcellCD4ex",
              "TcellCD8rest", "TcellCD8act", "TcellCD8ex")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type_Tcell, type, cell_type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type_Tcell %ni% c("other")) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type_Tcell = factor(cell_type_Tcell, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"high",
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))


p2 <- ggplot(frac_focus %>% filter(type=="IDH_O"),
             aes(x = cat_grade,
                 y = fraction,
                 fill = cell_type_Tcell)) +
  geom_col(position = "stack") +
  facet_grid(~ cell_type + type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_Tcell,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    y    = "Relative cell-type abundance",
    x    = "Cat_grade",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

##Significance CD4
keep_cd4 <- c("TcellCD4rest", "TcellCD4act", "TcellCD4ex")

# 1) per-sample fractions for CD4
frac_cd4 <- cells %>%
  filter(cell_type_Tcell %in% keep_cd4) %>%
  filter(type == "IDH_A") %>% 
  count(sample, section_grade, type, cell_type_Tcell) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# 2) pairwise Wilcoxon per CD4 subset, approximate p-values
stats_cd4 <- frac_cd4 %>%
  group_by(cell_type_Tcell) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(cell_type_Tcell, comparison)

stats_cd4 %>% arrange(p.value)

## Significance CD8
keep_cd8 <- c("TcellCD8rest", "TcellCD8act", "TcellCD8ex")

# 1) per-sample fractions for CD8 subsets
frac_cd8 <- cells %>%
  filter(cell_type_Tcell %in% keep_cd8) %>%
  filter(type == "IDH_A") %>% 
  count(sample, section_grade, type, cell_type_Tcell) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# 2) pairwise Wilcoxon per CD8 subset, approximate p-values
stats_cd8 <- frac_cd8 %>%
  group_by(cell_type_Tcell) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(cell_type_Tcell, comparison)

# View CD8 results
stats_cd8 %>% arrange(p.value)


# Vasc --------------------------------------------------------------------
keep_cts <- c("VascBBB", "Pericyte", "VascAng")

frac_per_grade <- cells %>%
  filter(!cell_type2 %in% c("low", "excluded")) %>% 
  count(section_grade, cell_type2, type) %>% 
  group_by(section_grade, type) %>% 
  mutate(fraction = n / sum(n)) %>% 
  ungroup()

frac_focus <- frac_per_grade %>%
  # keep only your 4 cell types
  filter(cell_type2 %in% keep_cts) %>%
  # set factor levels so stacking follows InfMg → GAM → MacScav → MacBorder
  mutate(cell_type2 = factor(cell_type2, levels = rev(keep_cts))) %>% 
  mutate(cat_grade = case_when(
    section_grade == 2 ~"low",
    section_grade == 3 ~"high",
  )) %>% 
  mutate(cat_grade = factor(cat_grade, levels = c("low","mid","high")))

p3 <- ggplot(frac_focus %>% filter(type=="IDH_O"),
             aes(x = cat_grade,
                 y = fraction,
                 fill = cell_type2)) +
  geom_col(position = "stack") +
  facet_grid(~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = an_col_vasc,
    breaks = rev(keep_cts)  # ensures legend matches stacking
  ) +
  labs(
    y    = "Relative cell-type abundance",
    x    = "Cat_grade",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.title        = element_text(size = 18),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background  = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text        = element_text(size = 13),
    legend.text       = element_text(size = 13)
  )

## Significance
keep_cts <- c("VascBBB", "Pericyte", "VascAng")

# 1) Per-sample fractions for your vascular subsets (IDH_A only)
frac_vasc <- cells %>%
  filter(cell_type2 %in% keep_cts,
         type == "IDH_A") %>%
  count(sample, section_grade, cell_type2) %>%      # count per sample×section_grade×cell_type2
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%              # fraction within each sample
  ungroup()

# 2) Pairwise Wilcoxon (2 vs 3, 3 vs 4, 2 vs 4) per vascular subset
stats_vasc <- frac_vasc %>%
  group_by(cell_type2) %>%
  do({
    d <- .
    r23 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,3)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_3")
    r34 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(3,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "3_vs_4")
    r24 <- broom::tidy(
      wilcox.test(frac ~ section_grade,
                  data  = filter(d, section_grade %in% c(2,4)),
                  exact = FALSE)
    ) %>% mutate(comparison = "2_vs_4")
    bind_rows(r23, r34, r24)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%  # BH-correction across all tests
  arrange(cell_type2, comparison)

# 3) View results
stats_vasc %>% arrange(p.value)
# Plot all ----------------------------------------------------------------
p1 <- p1 + labs(title = "Myeloid")
p2 <- p2 + labs(title = "T cell")
p3 <- p3 + labs(title = "Vasc")

(p1 + p2 + p3) + plot_layout(widths = c(1, 2, 1)) &
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

