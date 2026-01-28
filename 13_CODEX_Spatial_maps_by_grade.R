# Design for facet plots --------------------------------------------------
library(ggh4x)

design <- matrix(c(
  # Row 1 (low):  IDHA03 IDHA04 IDHA05 IDHA21 IDHO01 IDHO02 IDHO09 IDHA06 IDHO03 IDHO04 IDHO08 IDHO07 IDHO05
  3,    4,    5,   18, 6,   21,   22,   29,   23,   24,   28,   27,   25,
  # Row 2 (mid):  IDHA01 IDHA02 IDHA20 IDHA07 IDHA18 IDHA19 IDHA10 IDHA11 IDHA13 IDHA14 IDHA24 IDHA23
  1,    2,   17,    7,   15,   16,   10,   11,   13,   14,   20,   19,   NA,
  # Row 3 (high): IDHA08 IDHA09 IDHA12 IDHO10 IDHO06
  8,    9,   12,   30,   26,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA
), nrow = 3, byrow = TRUE)

# All cell trypes ---------------------------------------------------------
ggplot(cells %>% filter(cell_type %ni% c("low", "excluded")),
       aes(x = centroid_x, y = centroid_y, color = cell_type)) +
  geom_point(size = 0.01, alpha = 0.4) +
  scale_color_manual(values = an_col_simple) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  labs(color = "Cell type") +
  scale_y_reverse() +
  theme_minimal() +
  facet_manual(~type_sample, design = design)

p1 <- ggplot(
  cells %>% filter(cell_type %ni% c("low", "excluded")),
  aes(x = centroid_x / 1000, y = centroid_y / 1000, color = cell_type)  # µm -> mm
) +
  scattermore::geom_scattermore(pointsize = 1.5, alpha = 0.8) +
  # ggrastr::geom_point_rast(size = 0.001, alpha = 0.4, raster.dpi = 300) +
  scale_color_manual(values = an_col_simple) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  labs(
    color = "Cell Type",
    x = "Position x (mm)",
    y = "Position y (mm)"
  ) +
  scale_y_reverse() +
  theme_minimal() +
  ggh4x::facet_manual(~type_sample, design = design) +
  theme(
    panel.grid = element_blank(),                  
    panel.border = element_rect(                  
      color = "black", fill = NA, linewidth = 0.25
    ),
    legend.position = "bottom"
  )

# CairoPDF(
#   file = "cells_spatial_cairo.pdf",
#   width = 15,
#   height = 6, family = "Arial"
# )
# print(p1)
# dev.off()

# Spatial Map - IvyGAP ----------------------------------------------------
p2 <- ggplot(cells %>% filter(ivygap %ni% c("artifact")),
       aes(x = centroid_x / 1000, y = centroid_y / 1000, color = ivygap)  # µm -> mm
) +
  # ggrastr::geom_point_rast(size = 0.001, alpha = 0.4, raster.dpi = 300) +
  scattermore::geom_scattermore(pointsize = 1.5, alpha = 0.8) +
  scale_color_manual(values = an_col_ivygap) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  labs(
    color = "IvyGAP",
    x = "Position x (mm)",
    y = "Position y (mm)"
  ) +
  scale_y_reverse() +
  theme_minimal() +
  ggh4x::facet_manual(~type_sample, design = design) +
  theme(
    panel.grid = element_blank(),                  
    panel.border = element_rect(                  
      color = "black", fill = NA, linewidth = 0.25
    ),
    legend.position = "bottom"
  )

# CairoPDF(
#   file = "cells_spatial_cairo_ivygap.pdf",
#   width = 15,
#   height = 6, family = "Arial"
# )
# print(p2)
# dev.off()

