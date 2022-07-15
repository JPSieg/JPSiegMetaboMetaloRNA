
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(png)
library(ggtext)
library(rsvg)

####Figure_3A####

list.files("Figures/Figure_3_degredation")

bitmap <- rsvg_raw('Figures/Figure_3_degredation/Figure_3A_ILP_layout.svg', width = 600)
Figure_3A <- ggdraw() + draw_image(bitmap)

####Figure_3B####

list.files("Figures/Figure_3_degredation")

bitmap <- rsvg_raw('Figures/Figure_3_degredation/Figure_3B_Secondary_structure.svg', width = 600)
Figure_3B <- ggdraw() + draw_image(bitmap)

####Organize Gua data####

list.files("Figures/Figure_3_degredation")

vector.files = paste("Figures/Figure_3_degredation/Gua_data_files", list.files("Figures/Figure_3_degredation/Gua_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

####Make Figure 3C####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("25 mM Free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"),
                      labels = c("25 mM free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"))

unique(df$BP)

df$BP = factor(df$BP,
               levels = c("SS",
                          "NC",
                          "WC"))

table(df %>% filter(Condition == "2 mM free Mg") %>%
  select(BP))


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
  scale_color_manual(values = c( "red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  theme_classic()+
  ylab("Degradation (counts/hour)") +
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

####Figure 3D####


unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

Figure_3D = ggplot(df %>% filter(Reactivity > -200), aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition)) +
  facet_wrap(~BP, nrow = 1) +
  geom_line(mapping = aes(group = N), color = "dimgrey") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  #stat_compare_means(comparisons = comparisons, method = "t.test") +
  scale_color_manual(values = c("red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  theme_classic()+
  ggtitle("Guanine aptamer") +
  ylab("Degradation (counts/hour)") +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),axis.text.x = element_text(color = "Black", size = 14,
                                                                                             angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 10),
        axis.title.x = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.position = "none",
        plot.title = element_text(color = "Black", size = 16))

Figure_3D



####Organize CPEB3 data####

list.files("Figures/Figure_3_degredation")

vector.files = paste("Figures/Figure_3_degredation/CPEB3_data_files", list.files("Figures/Figure_3_degredation/CPEB3_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

df$Condition = factor(df$Condition,
                      levels = c("25 mM Free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"),
                      labels = c("25 mM free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"))

unique(df$BP)

df$BP = factor(df$BP,
               levels = c("SS",
                          "NC",
                          "WC"))

unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

df = df %>% filter(!is.na(Condition))

####Figure 3E####

table(df %>% filter(Condition == "2 mM free Mg") %>%
        select(BP))


Figure_3E = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  geom_line(mapping = aes(group = N), color = "dimgrey") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  #stat_compare_means(comparisons = comparisons, method = "t.test") +
  scale_color_manual(values = c("red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  theme_classic()+
  ggtitle("CPEB3 ribozyme") +
  ylab("Degradation (counts/hour)") +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14,
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

vector.files = paste("Figures/Figure_3_degredation/tRNA_data_files", list.files("Figures/Figure_3_degredation/tRNA_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

####Make Figure 3F####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("25 mM Free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"),
                      labels = c("25 mM free Mg",
                                 "2 mM free Mg",
                                 "Eco80",
                                 "NTPCM",
                                 "WMCM"))

df$BP = factor(df$BP,
               levels = c("SS",
                          "NC",
                          "WC"))

unique(df$Condition)

comparisons = list(c("25 mM free Mg", "2 mM free Mg"),
                   c("25 mM free Mg", "Eco80"),
                   c("25 mM free Mg", "NTPCM"),
                   c("25 mM free Mg", "WMCM"))

table(df %>% filter(Condition == "2 mM free Mg") %>%
        select(BP))

Figure_3F = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  geom_line(mapping = aes(group = N), color = "dimgrey") +
  facet_wrap(~BP) +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  #stat_compare_means(comparisons = comparisons, method = "t.test") +
  scale_color_manual(values = c("red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  theme_classic()+
  ggtitle("tRNAphe") +
  ylab("Degradation (counts/hour)") +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14,
                                   angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 10),
        axis.title.x = element_blank(),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),legend.text = element_text(color = "Black", size = 10),
        legend.position = "none",
        plot.title = element_text(color = "Black", size = 16))

####Check significance####

Figure_3D
Figure_3E
Figure_3F

####Consolidate plots into one plot####

Figure_3 = plot_grid(Figure_3A, Figure_3B, Figure_3C, Figure_3D, Figure_3E, Figure_3F,
                     labels = c("A", "B", "C", "D", "E", "F"),
                     ncol = 2, rel_heights = c(1, 1.4, 1.4),
                     label_size = 20)

#ggsave("Figures/Figure_3_degredation/Figure_3.svg", Figure_3,
#       width = 5, height = 6, scale = 2, bg = "white")
