
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(png)
library(ggtext)

####Organize Gua data####

vector.files = paste("Figures/Figure_2/Gua_data_files", list.files("Figures/Figure_2/Gua_data_files"), sep = "/")

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
  annotate("text", x = 32, y = -150, label = "P2", color = "white") +
  annotate("text", x = 47.5, y = -150, label = "P3", color = "white") +
  annotate("text", x = 60.5, y = -150, label = "P3", color = "white") +
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


unique(df$Condition)

comparisons = list(c("2 mM free Mg", "25 mM Free Mg"),
                   c("2 mM free Mg", "NTPCM"),
                   c("2 mM free Mg", "WMCM"),
                   c("2 mM free Mg", "Eco80"))

Figure_2D = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
  geom_boxplot(alpha = 0.01) +
  geom_beeswarm() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic()+
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
        legend.position = "none")



####Figure 2B####

list.files("Figures/Figure_2/2_structures")

Figure_2B = ggplot() +
  draw_image("Figures/Figure_2/2_structures/No_reactivity.png") +
  theme(panel.background = element_blank())

####Figure 2E####

list.files("Figures/Figure_2/2_structures")

Figure_2E = ggplot() +
  draw_image("Figures/Figure_2/2_structures/CPEB3_no_reactivity.png") +
  theme(panel.background = element_blank())

####Figure 2H####

list.files("Figures/Figure_2/2_structures")

Figure_2H = ggplot() +
  draw_image("Figures/Figure_2/2_structures/tRNA_no_reactivity.png") +
  theme(panel.background = element_blank())

####Figure_2A####

list.files("Figures/Figure_2")

PNG.2A = ggplot() +
  draw_image("Figures/Figure_2/Optimized_resolution.png") +
  theme(panel.background = element_blank())

####Organize CPEB3 data####

vector.files = paste("Figures/Figure_2/CPEB3_data_files", list.files("Figures/Figure_2/CPEB3_data_files"), sep = "/")

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

Figure_2F = ggplot(df, aes(x = N, y = Reactivity,
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


unique(df$Condition)

comparisons = list(c("2 mM free Mg", "25 mM Free Mg"),
                   c("2 mM free Mg", "NTPCM"),
                   c("2 mM free Mg", "WMCM"),
                   c("2 mM free Mg", "Eco80"))

Figure_2G = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
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
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),legend.text = element_text(color = "Black", size = 10),
        legend.position = "none")

####Organize tRNA data####

vector.files = paste("Figures/Figure_2/tRNA_data_files", list.files("Figures/Figure_2/tRNA_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)

df$BP = NA

df$BP[which(df$Dotbracket == ".")] = "Single stranded"
df$BP[which(df$Dotbracket != ".")] = "Base paired"

####Make Figure 2I####

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

df$Condition = factor(df$Condition,
                      levels = c("2 mM free Mg",
                                 "NTPCM",
                                 "WMCM",
                                 "Eco80",
                                 "25 mM Free Mg"))

Figure_2I = ggplot(df, aes(x = N, y = Reactivity,
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


unique(df$Condition)

comparisons = list(c("2 mM free Mg", "25 mM Free Mg"),
                   c("2 mM free Mg", "NTPCM"),
                   c("2 mM free Mg", "WMCM"),
                   c("2 mM free Mg", "Eco80"))

Figure_2J = ggplot(df, aes(x = Condition, y = Reactivity,
                           ymin = Reactivity - SE.k, ymax = Reactivity + SE.k,
                           color = Condition, group = Condition)) +
  facet_wrap(~BP) +
  stat_compare_means(comparisons = comparisons, size = 2.5,
                     label = "p.format", method = "t.test") +
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
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),legend.text = element_text(color = "Black", size = 10),
        legend.position = "none")


####Consolidate plots into one plot####

Figure_2BCD = plot_grid(Figure_2B, Figure_2C, Figure_2D, labels = c("B", "C", "D"),
                        label_size = 20, nrow = 1, rel_widths = c(1, 1.2, 1.2), hjust = c(-0.5, 0.75, 0.75))
Figure_2EFG = plot_grid(Figure_2E, Figure_2F, Figure_2G, labels = c("E", "F", "G"),
                        label_size = 20, nrow = 1, rel_widths = c(1, 1.2, 1.2), hjust = c(-0.5, 0.75, 0.75))

Figure_2HIJ = plot_grid(Figure_2H, Figure_2I, Figure_2J, labels = c("E", "F", "G"),
                        label_size = 20, nrow = 1, rel_widths = c(1, 1.2, 1.2), hjust = c(-0.5, 0.75, 0.75))


Figure_2ABCDEFGHIJ = plot_grid(PNG.2A, Figure_2BCD, Figure_2EFG, Figure_2HIJ, labels = c("A"), ncol = 1, rel_heights = c(2,1,1,1), label_size = 20)


ggsave("Figures/Figure_2/Figure_2.png", Figure_2ABCDEFGHIJ, width = 7, height = 12, scale = 2)
