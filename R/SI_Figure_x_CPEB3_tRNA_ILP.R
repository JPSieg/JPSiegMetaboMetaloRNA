
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(png)
library(ggtext)
library(rsvg)

#### SI Figure A CPEB3 secondary structure ####

list.files("Figures/SI_Figure_X_CPEB3_tRNAphe_ILP")

bitmap <- rsvg_raw('Figures/SI_Figure_X_CPEB3_tRNAphe_ILP/CPEB3_secondary_structure.svg', width = 600)
Figure_A <- ggdraw() + draw_image(bitmap)

Figure_A

#### SI Figure B Degradation verses N CPEB3 ####

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
               levels = c("Single stranded",
                          "Non-cannonical",
                          "WC",
                          "WC + Non-cannonical"))

df = df %>% filter(!is.na(Reactivity)) %>%
  filter(!is.na(Condition))

SI_figure_B = ggplot(df, aes(x = N, y = Reactivity,
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
  geom_segment(aes(x = 30, y = -150, xend = 36, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 36, y = -150, xend = 39, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 39, y = -150, xend = 45, yend = -150),
               color ="black", size = 5) +
  geom_segment(aes(x = 45, y = -150, xend = 50, yend = -150),
               color ="black", size = 0.5) +
  geom_segment(aes(x = 50, y = -150, xend = 56, yend = -150),
               color ="black", size = 5) +
  annotate("text", x = 28, y = -90, label = "P3-3'", color = "black") +
  annotate("text", x = 33, y = -150, label = "P1-3'", color = "white") +
  annotate("text", x = 42, y = -150, label = "P4-5'", color = "white") +
  annotate("text", x = 53, y = -150, label = "P4-3'", color = "white") +
  scale_color_manual(values = c("red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
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
        legend.text = element_text(color = "Black", size = 12),
        legend.position = c(0.5, 0.65),
        legend.title = element_blank(),
        legend.background = element_blank())

SI_figure_B

#### SI Figure A tRNA secondary structure ####

list.files("Figures/SI_Figure_X_CPEB3_tRNAphe_ILP")

bitmap <- rsvg_raw('Figures/SI_Figure_X_CPEB3_tRNAphe_ILP/tRNA_secondary_structure.svg', width = 600)
Figure_C <- ggdraw() + draw_image(bitmap)

Figure_C

#### SI Figure B Degradation verses N CPEB3 ####

vector.files = paste("Figures/Figure_3_degredation/tRNA_data_files", list.files("Figures/Figure_3_degredation/tRNA_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)

df = bind_rows(list.df)


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
               levels = c("Single stranded",
                          "Non-cannonical",
                          "WC",
                          "WC + Non-cannonical"))


SI_figure_D = ggplot(df, aes(x = N, y = Reactivity,
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
  scale_color_manual(values = c("red", "dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
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
        legend.text = element_text(color = "Black", size = 12),
        legend.position = c(0.5, 0.65),
        legend.title = element_blank(),
        legend.background = element_blank())

SI_figure_D


SI_figure_X = plot_grid(Figure_A, SI_figure_B, Figure_C, SI_figure_D, labels = c("A", "B", "C", "D"),
          label_size = 20)

list.files("Figures/SI_Figure_X_CPEB3_tRNAphe_ILP")

ggsave("Figures/SI_Figure_X_CPEB3_tRNAphe_ILP/SI_figure_x_CPEB3_tRNA_ILP.svg",
       SI_figure_X,
       bg = "white",
       scale = 3)
ggsave("Figures/SI_Figure_X_CPEB3_tRNAphe_ILP/SI_figure_x_CPEB3_tRNA_ILP.png",
       SI_figure_X,
       bg = "white",
       scale = 3)
