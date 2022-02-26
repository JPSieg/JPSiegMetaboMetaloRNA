library(tidyverse)
library(cowplot)

PNG.2mMFree = ggplot() +
  draw_image("Figures/Figure_2/2_structures/2mMFree.png") +
  theme(panel.background = element_blank())
PNG.NTPCM = ggplot() +
  draw_image("Figures/Figure_2/2_structures/NTPCM.png") +
  theme(panel.background = element_blank())
PNG.WMCM = ggplot() +
  draw_image("Figures/Figure_2/2_structures/WMCM.png") +
  theme(panel.background = element_blank())
PNG.Eco80 = ggplot() +
  draw_image("Figures/Figure_2/2_structures/Eco80.png") +
  theme(panel.background = element_blank())
PNG.25mMFree = ggplot() +
  draw_image("Figures/Figure_2/2_structures/25mMFree.png") +
  theme(panel.background = element_blank())
PNG.Legend = ggplot() +
  draw_image("Figures/Figure_2/2_structures/Legend.png") +
  theme(panel.background = element_blank())

Figure_2B = plot_grid(PNG.2mMFree, PNG.NTPCM, PNG.WMCM,
                      PNG.Eco80, PNG.25mMFree, PNG.Legend,
                      labels = c("2 mM Free", "NTPCM", "WMCM", "Eco80", "25 mM free", ""))

ggsave("Figures/SI_figure_x_reactivity_mapped/Reactivity_mapped.png",
       Figure_2B)
