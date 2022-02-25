
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(png)

####Organize data####

vector.files = paste("Figures/Figure_2/Data_files", list.files("Figures/Figure_2/Data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

df$BP = NA

df$BP[which(df$Dotbracket == ".")] = "Single stranded"
df$BP[which(df$Dotbracket != ".")] = "Base paired"

####Make Figure 2C####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("2 mM free Mg",
                                 "NTPCM",
                                 "WMCM",
                                 "Eco80",
                                 "25 mM Free Mg"))

Figure_2C = ggplot(df, aes(x = N, y = Reactivity,
               ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
               color = Condition, group = Condition)) +
  geom_pointrange() +
  geom_line() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ylab("Estimated dCounts/dt") +
  xlab("Nucleotide") +
  xlim(29, 63) +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 8),
        legend.position = c(0.55, 0.85),
        legend.title = element_blank(),
        legend.background = element_blank())


unique(df$Condition)

comparisons = list(c("2 mM free Mg", "25 mM Free Mg"),
                   c("2 mM free Mg", "NTPCM"),
                   c("2 mM free Mg", "WMCM"),
                   c("2 mM free Mg", "Eco80"))

Figure_2D = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5) +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ylab("Estimated dCounts/dt") +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16,
                                   angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 10),
        axis.title.x = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.position = "none")



####Figure 2B####

list.files("Figures/Figure_2/2_structures")

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



####Figure_2A####

list.files("Figures/Figure_2")

PNG.2A = ggplot() +
  draw_image("Figures/Figure_2/Optimized_resolution.png") +
  theme(panel.background = element_blank())

Figure_2CD = plot_grid(Figure_2C, Figure_2D, ncol = 1, labels = c("C", "D"), label_x = -0.05, rel_heights = c(1, 1.3), label_size = 20)
Figure_2BCD = plot_grid(Figure_2B, Figure_2CD, rel_widths = c(1.5, 1), labels = c("B", ""), label_size = 20)
Figure_2ABCD = plot_grid(PNG.2A, Figure_2BCD, labels = c("A"), ncol = 1, rel_heights = c(5,4.5), label_size = 20)

ggsave("Figures/Figure_2/Figure_2.png", Figure_2ABCD, width = 7, height = 8, scale = 2)

