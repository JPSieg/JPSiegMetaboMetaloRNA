library(tidyverse)
library(cowplot)
library(ggtext)

list.files("Figures/Figure_4/")

####Experiment design image####

PNG.design = ggplot() +
  draw_image("Figures/Figure_4/Experimental_design.png") +
  theme(panel.background = element_blank())

####Rotate symbol####

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


####Figure 4B Rg and Dmax####

df.Rg = read.csv("Figures/Figure_4/Rg_data.csv")

df.Rg$Parameter = factor(df.Rg$Parameter,
                         levels = c("Gunier Rg", "DENSS Rg",
                                    "Pr Rg", "Dmax"),
                         labels = c("Gunier R<sub>g</sub>", "DENSS R<sub>g</sub>",
                                    "p(r) R<sub>g</sub>", "D<sub>max</sub>"))

df.Rg$Condition = factor(df.Rg$Condition,
                          levels = c("2 mM Mg",
                                     "NTPCM",
                                     "WMCM",
                                     "Eco80",
                                     "25 mM Mg"))

Figure_4B = ggplot(df.Rg, aes(x = Condition,
                  y = Value,
                  ymin = Value - Error,
                  ymax = Value + Error,
                  fill = Condition,
                  label = round(Value, digits = 1))) +
  scale_fill_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  geom_text(nudge_y = 6) +
  facet_wrap(~Parameter) +
  ylab("Size (\u212b)") +
  coord_cartesian(ylim = c(20, 90)) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.background = element_blank())


####Consolidate plots####

Figure_4AB = plot_grid(PNG.design, Figure_4B, ncol = 1,
                       rel_heights = c(2, 1), labels = c("A", "B"),
                       label_size = 20)

Figure_4C = plot_grid(F2mMMg,
                     FNTPCM,
                     FWMCM,
                     FEco80,
                     F2mMMg,
                     ncol = 1,
                     labels = c("2 mM Free",
                                "2 mM Free + NTPCM",
                                "2 mM Free + WMCM",
                                "2 mM Free + Eco80",
                                "25 mM Free"),
                     label_size = 18,
                     label_x = c(0.25, 0.1, 0.1, 0.1, 0.25))

Figure_4ABC = plot_grid(Figure_4AB, Figure_4C,
                        nrow = 1, labels = c("", "C"),
                        label_size = 20)

ggsave("Figures/Figure_4/Figure_4.png", Figure_4ABC, scale = 2)
