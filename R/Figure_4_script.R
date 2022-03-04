library(tidyverse)
library(cowplot)
library(PNG)

list.files("Figures/Figure_4/Images_from_pymol")



PNG.rotate = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/Rotate_symbol.png") +
  theme(panel.background = element_blank())


####2mM Mg####

PNG.2mM0 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/2mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())
PNG.2mM90 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/2mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

F2mMMg = plot_grid(PNG.2mM0, PNG.2mM90, nrow = 1) +
  draw_image("Figures/Figure_4/Images_from_pymol/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

####NTPCM####

PNG.NTPCM0 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/NTPCM_angle_0deg.png") +
  theme(panel.background = element_blank())
PNG.NTPCM90 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/NTPCM_angle_90deg.png") +
  theme(panel.background = element_blank())

FNTPCM =plot_grid(PNG.NTPCM0, PNG.NTPCM90, nrow = 1)+
  draw_image("Figures/Figure_4/Images_from_pymol/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

####WMCM####

PNG.WMCM0 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/WMCM_angle_0deg.png") +
  theme(panel.background = element_blank())
PNG.WMCM90 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/WMCM_angle_90deg.png") +
  theme(panel.background = element_blank())

FWMCM =plot_grid(PNG.WMCM0, PNG.WMCM90, nrow = 1) +
  draw_image("Figures/Figure_4/Images_from_pymol/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())


####Eco80####

PNG.Eco800 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/Eco80_angle_0deg.png") +
  theme(panel.background = element_blank())
PNG.Eco8090 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/Eco80_angle_90deg.png") +
  theme(panel.background = element_blank())

FEco80 =plot_grid(PNG.Eco800, PNG.Eco8090, nrow = 1) +
  draw_image("Figures/Figure_4/Images_from_pymol/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

####25mM Mg####

PNG.25mM0 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/25mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())
PNG.25mM90 = ggplot() +
  draw_image("Figures/Figure_4/Images_from_pymol/25mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

F25mMMg =plot_grid(PNG.25mM0, PNG.25mM90, nrow = 1)

####Consolidate plots####

Figure_4 = plot_grid(F2mMMg,
                     FNTPCM,
                     FWMCM,
                     FEco80,
                     F2mMMg,
                     ncol = 1,
                     labels = c("2 mM Free",
                                "2 mM Free + WMCM",
                                "2 mM Free + NTPCM",
                                "2 mM Free + Eco80",
                                "25 mM Free"), label_size = 5)

ggsave("Figures/Figure_4/Figure_4.png", width = 2)
