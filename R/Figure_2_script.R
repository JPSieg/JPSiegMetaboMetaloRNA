
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(png)
library(ggtext)
library(rsvg)

####Organize Gua data####

vector.files = paste("Figures/Figure_2/Gua_data_files", list.files("Figures/Figure_2/Gua_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

####Make Figure 2C####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM Free Mg"),
                      labels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM free Mg"))

unique(df$BP)

df$BP = factor(df$BP,
               levels = c("Single stranded",
                          "Non-cannonical",
                          "WC",
                          "WC + Non-cannonical"))


Figure_3C = ggplot(df, aes(x = N, y = Reactivity,
               ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
               color = Condition, group = Condition)) +
  geom_pointrange() +
  geom_line() +
  geom_segment(aes(x = 29, y = -150, xend = 30, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 30, y = -150, xend = 34, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 34, y = -150, xend = 45, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 45, y = -150, xend = 50, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 50, y = -150, xend = 58, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 58, y = -150, xend = 63, yend = -150),
               color ="black", size = 5) +
  annotate("text", x = 32, y = -150, label = "P2-3'", color = "white") +
  annotate("text", x = 47.5, y = -150, label = "P3-5'", color = "white") +
  annotate("text", x = 60.5, y = -150, label = "P3-3'", color = "white") +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ylab("Estimated dCounts/dt") +
  xlab("Nucleotide") +
  xlim(29, 63) +
  ylim(-200, 1500) +
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

Figure_3C


unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

Figure_3D = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP, nrow = 1) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ggtitle("Guanine aptamer") +
  ylab("Estimated dCounts/dt") +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),axis.text.x = element_text(color = "Black", size = 16,
                                   angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 10),
        axis.title.x = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.position = "none",
        plot.title = element_text(color = "Black", size = 16))

Figure_3D

####Figure_3A####

list.files("Figures/Figure_2")

bitmap <- rsvg_raw('Figures/Figure_2/Figure_3A_ILP_layout.svg', width = 600)
Figure_3A <- ggdraw() + draw_image(bitmap)

####Figure_3B####

list.files("Figures/Figure_2")

bitmap <- rsvg_raw('Figures/Figure_2/Figure_3B_Secondary_structure.svg', width = 600)
Figure_3B <- ggdraw() + draw_image(bitmap)

####Organize CPEB3 data####

vector.files = paste("Figures/Figure_2/CPEB3_data_files", list.files("Figures/Figure_2/CPEB3_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

df$Condition = factor(df$Condition,
                      levels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM Free Mg"),
                      labels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM free Mg"))

unique(df$BP)

df$BP = factor(df$BP,
               levels = c("Single stranded",
                          "Non-cannonical",
                          "WC",
                          "WC + Non-cannonical"))

####Make SI Figure for CPEB3####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

SI_figure_X.1 = ggplot(df, aes(x = N, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  geom_pointrange() +
  geom_line() +
  geom_segment(aes(x = 22, y = -150, xend = 27, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 27, y = -150, xend = 29, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 29, y = -150, xend = 31, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 31, y = -150, xend = 36, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 36, y = -150, xend = 39, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 39, y = -150, xend = 45, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 45, y = -150, xend = 50, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 50, y = -150, xend = 56, yend = -150),
               color ="black", size = 5) +
  annotate("text", x = 28, y = -150, label = "P3", color = "white") +
  annotate("text", x = 33.5, y = -150, label = "P1", color = "white") +
  annotate("text", x = 41.5, y = -150, label = "P4", color = "white") +
  annotate("text", x = 53, y = -150, label = "P4", color = "white") +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ylab("Estimated dCounts/dt") +
  xlab("Nucleotide") +
  xlim(22, 56) +
  ylim(-200, 1500) +
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

SI_figure_X.1

unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

Figure_3E = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ggtitle("CPEB3 ribozyme") +
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
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),legend.text = element_text(color = "Black", size = 10),
        legend.position = "none",
        plot.title = element_text(color = "Black", size = 16))

Figure_3E

####Organize tRNA data####

vector.files = paste("Figures/Figure_2/tRNA_data_files", list.files("Figures/Figure_2/tRNA_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

####Make Figure 2I####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM Free Mg"),
                      labels = c("2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM",
                                 "25 mM free Mg"))

df$BP = factor(df$BP,
               levels = c("Single stranded",
                          "Non-cannonical",
                          "WC",
                          "WC + Non-cannonical"))


SI_figure_X.2 = ggplot(df, aes(x = N, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  geom_pointrange() +
  geom_line() +
  geom_segment(aes(x = 24, y = -300, xend = 25, yend = -300),
               color ="black", size = 5) +
  geom_segment(aes(x = 25, y = -300, xend = 27, yend = -300),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 27, y = -300, xend = 31, yend = -300),
               color ="black", size = 5) +
  geom_segment(aes(x = 31, y = -300, xend = 39, yend = -300),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 39, y = -300, xend = 43, yend = -300),
               color ="black", size = 5) +
  geom_segment(aes(x = 43, y = -300, xend = 49, yend = -300),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 49, y = -300, xend = 53, yend = -300),
               color ="black", size = 5) +
  geom_segment(aes(x = 53, y = -300, xend = 61, yend = -300),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 61, y = -300, xend = 65, yend = -300),
               color ="black", size = 5) +
  annotate("text", x = 29, y = -300, label = "P3-5'", color = "white") +
  annotate("text", x = 41, y = -300, label = "P3-3'", color = "white") +
  annotate("text", x = 51, y = -300, label = "P4-5'", color = "white") +
  annotate("text", x = 63, y = -300, label = "P4-3'", color = "white") +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ylab("Estimated dCounts/dt") +
  xlab("Nucleotide") +
  xlim(24, 65) +
  ylim(-400, 5000) +
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

SI_figure_X.2

unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

Figure_3F = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
  ggtitle("tRNA") +
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
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),legend.text = element_text(color = "Black", size = 10),
        legend.position = "none",
        plot.title = element_text(color = "Black", size = 16))

Figure_3F

####Consolidate plots into one plot####

Figure_3ABC = plot_grid(Figure_3A, Figure_3B, Figure_3C, nrow = 1, labels = c("A", "B", "C"), label_size = 20)

Figure_3DEF = plot_grid(Figure_3D, Figure_3E, Figure_3F, nrow = 1, labels = c("D", "E", "F"), rel_widths = c(1.7, 1, 1), label_size = 20)

Figure_3 = plot_grid(Figure_3ABC, Figure_3DEF, ncol = 1)

Figure_3

ggsave("Figures/Figure_2/Figure_3.svg", Figure_3, width = 7, height = 4, scale = 3.2, bg = "white")
