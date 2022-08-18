library(tidyverse)
library(MeltR)
library(viridis)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(ggpubr)


####Raw data plot####

Figure_2/Fluorescence_data

list.files("Figures/Figure_2_Table_3_SI_table_5_helix_stability")

df = read.csv("Figures/Figure_2_Table_3_SI_table_5_helix_stability/Fluorescence_data/js5060_Helix_F_formatted.csv")
df.Ks = read.csv("Figures/Figure_2_Table_3_SI_table_5_helix_stability/Fits_summary_Ks.csv")


head(df)
head(df.Ks)


df.Ks = df.Ks %>%
  filter(Condition == "Monovalent") %>%
  filter(Helix == "F")


readings = c(20, 40, 60, 67, 75, 80, 85)

Mmodel1 = function(x){df.Ks$Fmax[readings[1]] + (df.Ks$Fmin[readings[1]] - df.Ks$Fmax[readings[1]])*(((10^9)/df.Ks$K[readings[1]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[1]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel2 = function(x){df.Ks$Fmax[readings[2]] + (df.Ks$Fmin[readings[2]] - df.Ks$Fmax[readings[2]])*(((10^9)/df.Ks$K[readings[2]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[2]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel3 = function(x){df.Ks$Fmax[readings[3]] + (df.Ks$Fmin[readings[3]] - df.Ks$Fmax[readings[3]])*(((10^9)/df.Ks$K[readings[3]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[3]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel4 = function(x){df.Ks$Fmax[readings[4]] + (df.Ks$Fmin[readings[4]] - df.Ks$Fmax[readings[4]])*(((10^9)/df.Ks$K[readings[4]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[4]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel5 = function(x){df.Ks$Fmax[readings[5]] + (df.Ks$Fmin[readings[5]] - df.Ks$Fmax[readings[5]])*(((10^9)/df.Ks$K[readings[5]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[5]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel6 = function(x){df.Ks$Fmax[readings[6]] + (df.Ks$Fmin[readings[6]] - df.Ks$Fmax[readings[6]])*(((10^9)/df.Ks$K[readings[6]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[6]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}
Mmodel7 = function(x){df.Ks$Fmax[readings[7]] + (df.Ks$Fmin[readings[7]] - df.Ks$Fmax[readings[7]])*(((10^9)/df.Ks$K[readings[7]]+(200/1.3470043)+x)-((((10^9)/df.Ks$K[readings[7]]+(200/1.3470043)+x)^2)-(4*(200/1.3470043)*x))^(1/2))/(2*(200/1.3470043))}

df.raw = df %>% filter(Reading %in% readings)

unique(df.raw$Temperature)[1:13]

raw.plot = ggplot(df %>% filter(Reading %in% readings)) +
  geom_point(mapping = aes(x = B,
                           y = Emission,
                           color = factor(floor(Temperature)))) +
  stat_function(fun = Mmodel1, color = viridis(7, option = "H")[1]) +
  stat_function(fun = Mmodel2, color = viridis(7, option = "H")[2]) +
  stat_function(fun = Mmodel3, color = viridis(7, option = "H")[3]) +
  stat_function(fun = Mmodel4, color = viridis(7, option = "H")[4]) +
  stat_function(fun = Mmodel5, color = viridis(7, option = "H")[5]) +
  stat_function(fun = Mmodel6, color = viridis(7, option = "H")[6]) +
  stat_function(fun = Mmodel7, color = viridis(7, option = "H")[7]) +
  scale_color_manual(values = viridis(7, option = "H")) +
  theme_classic() +
  ylab("FAM-RNA emission") +
  xlab("[RNA-BHQ1] (nM)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = "none")

####Make vh plots####

df.Ks = read.csv("Figures/Figure_2_Table_3_SI_table_5_helix_stability/Fits_summary_Ks.csv")

df.Ks$Condition = factor(df.Ks$Condition,
                         levels = c("Monovalent","Ecoli80", "NTPCM", "WMCM"),
                         labels = c("Monovalent","Eco80", "NTPCM", "WMCM"))


df.Ks = df.Ks %>% filter(Helix == "F")


vh.plot = ggplot(df.Ks %>%
                    filter(In_Kd_range == TRUE) %>%
                    filter(SE.lnK <= K_error),
                  aes(x = invT, y = lnK,
                      ymin = lnK - SE.lnK,
                      ymax = lnK +SE.lnK,
                      color = Condition,
                      shape = Condition)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(6, 3, 1)])) +
  geom_pointrange() +
  ylab("ln[K (1/M)]") +
  xlab("1/Temperature (K)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.3),
        legend.background = element_blank())

####Make thermo plot####

df.vh = read.csv("Figures/Figure_2_Table_3_SI_table_5_helix_stability/Fits_summary_vh.csv")
head(df.vh)


df.vh$Condition = factor(df.vh$Condition,
                         levels = c("Monovalent","Ecoli80", "NTPCM", "WMCM"),
                         labels = c("Monovalent","Eco80", "NTPCM", "WMCM"))
df.vh$Helix = factor(df.vh$Helix,
                     levels = c("I", "F", "J", "G", "H"),
                     labels = c("1:CGGAUGGC/GCCAUCCG",
                                "2:CGCAUCCU/AGGAUGCG",
                                "3:CGUAUGUA/UACAUACG",
                                "4:CCAUAUCA/UGAUAUGG",
                                "5:CCAUAUUA/UAAUAUGG"))


df.vh$K = exp(-as.numeric(df.vh$G)/(0.00198720425864083*(273.15 + 37)))
df.vh$SE.K = df.vh$K*as.numeric(df.vh$SE.G)/as.numeric(df.vh$G)
dG.plot = ggplot(data = df.vh %>% filter(Method == "1 VH plot"),
       mapping = aes(x = Condition, y = G, ymin = G - 0.015*G, ymax = G + 0.015*G)) +
  facet_wrap(~Helix, nrow = 1, scales = "free") +
  geom_point(stat="identity", size = 5) +
  geom_errorbar(size = 1.5) +
  theme_classic() +
  ylab("K at 37\u00b0C (1/M)") +
  scale_y_reverse(breaks = scales::pretty_breaks(n = 5)) +
  #coord_cartesian(ylim=c(5,12.5)) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        legend.title = element_text(color = "Black", size = 14))

dG.plot

####Read in figure 3A####

library(svglite)
library(magick)
library(rsvg)
library(grobblR)
library(grid)

list.files("Figures/Figure_2_Table_3_SI_table_5_helix_stability")

bitmap <- rsvg_raw('Figures/Figure_2_Table_3_SI_table_5_helix_stability/Figure_3A.svg', width = 600)
image <- ggdraw() + draw_image(bitmap)

####Make figure 3####

Figure_3ABC = plot_grid(image, raw.plot, vh.plot, nrow = 1,
                        labels = c("A", "B", "C"))

Figure_ABCDE = plot_grid(Figure_3ABC,
                         dG.plot,
                         ncol = 1,
                         labels = c("", "D"))

Figure_ABCDE

ggsave("Figures/Figure_2_Table_3_SI_table_5_helix_stability/Figure_2.svg",
       Figure_ABCDE,
       scale = 2.5,
       width = 7, height = 4.5, units = "in",
       bg = "white")
