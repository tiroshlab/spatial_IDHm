# Levels for Compositions -------------------------------------------------
levels_cell_type <- c(
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
  "cAC2",
  "cAC1",
  "cOPC",
  "cGEM"
)

order_legend_cell_type <-  c(
  "cMESHyp",
  "cMES",
  "cUndiff",
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

malig <- c("cMESHyp", 
           "cMES", 
           "cAC1",
           "cAC2",
           "cOPC",
           "cGEM",
           "cUndiff")

gcells <- gcells %>% 
  mutate(grade_class = case_when(
    sample %in% samples_low ~ "low",
    sample %in% samples_mid ~ "mid",
    sample %in% samples_high ~ "high"
    ))

# Composition bar plot per type ------------------------------------------------------
## Per type and optionally only for "ct"
frac_per_type <- gcells %>%
  filter(cell_type %ni% c("low", "excluded")) %>% 
  filter(ivygap == "ct") %>%
  count(type, cell_type) %>%      
  group_by(type) %>%               
  mutate(fraction = n / sum(n)) %>%  
  ungroup()

frac_per_type <- frac_per_type %>% 
  mutate(fraction_neg = if_else(cell_type %in% malig,
                                true = fraction,
                                false = fraction*-1)) %>% 
  group_by(type) %>% 
  mutate(order = sum(fraction_neg)) %>% 
  ungroup()

ggplot(frac_per_type, aes(x=type, 
                            y=fraction_neg, 
                            fill=factor(cell_type, levels = levels_cell_type
                            ))) +
  geom_bar(stat="identity") +  
  theme_minimal() + 
  theme(axis.text.x = element_text(size= 12, angle=45, hjust = 1), 
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  scale_fill_manual(values = an_col_simple, breaks = order_legend_cell_type)+
  theme(legend.text=element_text(size=13)) +
  geom_hline(yintercept = 0) +
  labs(y="relative cell type abundance", x="type", fill="cell type")

##Significance:
frac_sample <- gcells %>%
  filter(!(cell_type %in% c("low", "excluded"))) %>%
  filter(ivygap == "ct") %>%
  count(sample, type, cell_type) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

frac_sample %>%
  group_by(cell_type) %>%
  do(broom::tidy(wilcox.test(frac ~ type, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj)

# Composition bar plot per sample -------------------------------------------
frac_per_sample <- gcells %>%
  filter(cell_type %ni% c("low", "excluded")) %>% 
  count(sample, cell_type, type) %>%       # Count gcells for each sample and cell type
  group_by(sample) %>%               # Group by sample
  mutate(fraction = n / sum(n)) %>%  # Calculate fraction for each cell type
  ungroup()                        # Ungroup the data

frac_per_sample <- frac_per_sample %>% 
  mutate(fraction_neg = if_else(cell_type %in% malig,
                                true = fraction,
                                false = fraction*-1)) %>% 
  group_by(sample, type) %>% 
  mutate(order = sum(fraction_neg)) %>% 
  ungroup()

ggplot(frac_per_sample, aes(
  x    = fct_reorder(sample, order),
  y    = fraction_neg,
  fill = factor(cell_type, levels = levels_cell_type)
)) +
  geom_col() +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = an_col_simple, breaks = order_legend_cell_type) +
  geom_hline(yintercept = 0) +
  labs(
    y    = "relative cell type abundance",
    x    = "sample",
    fill = "cell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x         = element_text(size = 12, angle = 45, hjust = 1),
    axis.title           = element_text(size = 18),
    plot.title           = element_text(size = 18),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    ## draw a box around each facet:
    panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.spacing        = unit(0.5, "lines"),
    ## (optional) box in facet strips too:
    strip.background     = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text           = element_text(size = 13),
    legend.text          = element_text(size = 13)
  )

# Composition per grade and type --------------------------------------------
frac_per_grade <- gcells %>%
  filter(cell_type %ni% c("low", "excluded")) %>% 
  count(section_grade, cell_type, type) %>%       # Count gcells for each sample and cell type
  group_by(section_grade, type) %>%               # Group by sample
  mutate(fraction = n / sum(n)) %>%  # Calculate fraction for each cell type
  ungroup()                        # Ungroup the data

frac_per_grade <- frac_per_grade %>% 
  mutate(fraction_neg = if_else(cell_type %in% malig,
                                true = fraction,
                                false = fraction*-1)) %>% 
  group_by(section_grade, type) %>% 
  mutate(order = sum(fraction_neg)) %>% 
  ungroup()

ggplot(frac_per_grade, aes(
  x    = factor(section_grade, levels = c(2:4)),
  y    = fraction_neg,
  fill = factor(cell_type, levels = levels_cell_type)
)) +
  geom_col() +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = an_col_simple, breaks = order_legend_cell_type) +
  geom_hline(yintercept = 0) +
  labs(
    y    = "relative cell type abundance",
    x    = "Grade",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x         = element_text(size = 12, angle = 45, hjust = 1),
    axis.title           = element_text(size = 18),
    plot.title           = element_text(size = 18),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    ## draw a box around each facet:
    panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.spacing        = unit(0.5, "lines"),
    ## (optional) box in facet strips too:
    strip.background     = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text           = element_text(size = 13),
    legend.text          = element_text(size = 13)
  )

## Significance IDH-A only
# 1) Compute per‐sample fractions (just like before)
frac_sample <- gcells %>%
  filter(!(cell_type %in% c("low", "excluded", "Bcell", "Neutro", "TcellCD4", "TcellCD8"))) %>%
  filter(type == "IDH_A") %>% 
  count(sample, section_grade, type, cell_type) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# 2) Pairwise Wilcoxon for section_grades 2→3 and 3→4, per cell_type
pairwise_section_grade_stats <- frac_sample %>%
  group_by(cell_type) %>%
  do(
    bind_rows(
      # 2 vs 3
      broom::tidy(
        wilcox.test(frac ~ section_grade,
                    data = filter(., section_grade %in% c(2,3)))
      ) %>% mutate(comparison = "2_vs_3"),
      # 3 vs 4
      broom::tidy(
        wilcox.test(frac ~ section_grade,
                    data = filter(., section_grade %in% c(3,4)))
      ) %>% mutate(comparison = "3_vs_4"),
      # 2 vs 4
      broom::tidy(
        wilcox.test(frac ~ section_grade,
                    data = filter(., section_grade %in% c(2,4)))
      ) %>% mutate(comparison = "2_vs_4"),
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p_adj)

# Inspect
pairwise_section_grade_stats %>% print(n=Inf)

# Composition bar plot per grade class --------------------------------------------
frac_per_grade <- gcells %>%
  filter(cell_type %ni% c("low", "excluded")) %>% 
  count(grade_class, cell_type) %>%       # Count gcells for each sample and cell type
  group_by(grade_class) %>%               # Group by sample
  mutate(fraction = n / sum(n)) %>%  # Calculate fraction for each cell type
  ungroup()                        # Ungroup the data

frac_per_grade <- frac_per_grade %>% 
  mutate(fraction_neg = if_else(cell_type %in% malig,
                                true = fraction,
                                false = fraction*-1)) %>% 
  group_by(grade_class) %>% 
  mutate(order = sum(fraction_neg)) %>% 
  ungroup()

ggplot(frac_per_grade, aes(
  x    = factor(grade_class, levels = c("low", "mid", "high")),
  y    = fraction_neg,
  fill = factor(cell_type, levels = levels_cell_type)
)) +
  geom_col() +
  # facet_grid(~grade_class, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = an_col_simple, breaks = order_legend_cell_type) +
  geom_hline(yintercept = 0) +
  labs(
    y    = "relative cell type abundance",
    x    = "Grade Class",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x         = element_text(size = 12, angle = 45, hjust = 1),
    axis.title           = element_text(size = 18),
    plot.title           = element_text(size = 18),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    ## draw a box around each facet:
    panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.spacing        = unit(0.5, "lines"),
    ## (optional) box in facet strips too:
    strip.background     = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text           = element_text(size = 13),
    legend.text          = element_text(size = 13)
  )

# Composition per IvyGAP ----------------------------------------------------
##Per IVYgap
frac_per_ivygap <- gcells %>%
  filter(cell_type %ni% c("low", "excluded")) %>% 
  filter(ivygap %ni% c("none", "artifact")) %>%
  count(ivygap, cell_type, type) %>%       # Count gcells for each ivygap and cell type
  group_by(ivygap, type) %>%               # Group by ivygap
  mutate(fraction = n / sum(n)) %>%  # Calculate fraction for each cell type
  ungroup()                        # Ungroup the data

frac_per_ivygap <- frac_per_ivygap %>% 
  mutate(fraction_neg = if_else(cell_type %in% malig,
                                true = fraction,
                                false = fraction*-1)) %>% 
  group_by(ivygap, type) %>% 
  mutate(order = sum(fraction_neg)) %>% 
  ungroup()

ggplot(frac_per_ivygap, aes(
  x    = fct_reorder(ivygap, order),
  y    = fraction_neg,
  fill = factor(cell_type, levels = levels_cell_type)
)) +
  geom_col(width = 0.8) +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = an_col_simple, breaks = order_legend_cell_type) +
  geom_hline(yintercept = 0) +
  labs(
    y     = "relative cell type abundance",
    x     = "IvyGAP",
    fill  = "cell type"
  ) +
  theme_minimal() +
  theme(
    # axes & text
    axis.text.x         = element_text(size = 12, angle = 45, hjust = 1),
    axis.title          = element_text(size = 18),
    plot.title          = element_text(size = 18, face = "bold", hjust = 0.5),
    
    # grid and background
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    
    # facet panels
    panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.spacing       = unit(0.5, "lines"),
    
    # facet strip (the grey titles)
    strip.background    = element_rect(colour = "black", fill = "grey90", linewidth = 0.7),
    strip.text          = element_text(size = 13),
    
    # legend
    legend.text         = element_text(size = 13)
  )

