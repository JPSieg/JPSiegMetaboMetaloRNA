library(tidyverse)
library(MeltR)
library(viridis)
library(ggbeeswarm)
library(cowplot)


####Raw data plot####

list.files("Figures/Figure_3/Fluorescence_data")

df = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_F_formatted.csv")
df.Ks = read.csv("Figures/Figure_3/Fits_summary_Ks.csv")


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
  ylab("Emission") +
  xlab("[BHQ1] (nM)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = "none")

raw.plot

#####Make Tms plot####

list.files("Figures/Figure_3")

df.Tms = read.csv("Figures/Figure_3/Fits_summary_Tms.csv")

df.Tms$Condition = factor(df.Tms$Condition,
                          levels = c("Monovalent", "NTPCM", "WMCM", "Ecoli80"))

head(df.Tms)

df.Tms$logB = log10(df.Tms$B)

df.Tms = df.Tms %>% filter(B >= 250)

ggplot(df.Tms, aes(x = logB, y = Tm, color = Condition)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Helix, scales = "free")

lm.fit = lm(Tm ~ logB * Helix * Condition, df.Tms)

df.Tms$predict = predict(lm.fit)

Tm.plot = ggplot() +
  scale_color_viridis(name = "[BHQ1] (nM)") +
  facet_wrap(~Helix, nrow = 1) +
  geom_beeswarm(data = df.Tms,
                mapping = aes(x = Condition, y = Tm, color = B)) +
  geom_line(data = df.Tms,
            mapping = aes(x = Condition, y = predict, color = B, group = logB)) +
  theme_classic() +
  ylab("Tm (\u00b0C)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = c(0.5, 0.6))

Tm.plot

####Make vh plots####

df.Ks = read.csv("Figures/Figure_3/Fits_summary_Ks.csv")

list.files("Figures/Figure_3")

head(df.Ks)

unique(df.Ks$Condition)

df.Ks$lnK = log(df.Ks$K)
df.Ks$invT = 1/(273.15 + df.Ks$Temperature)
df.Ks$SE.lnK = df.Ks$SE.K/df.Ks$K

df.Ks$Condition = factor(df.Ks$Condition,
                         levels = c("Monovalent", "NTPCM", "WMCM", "Ecoli80"))


df.Ks = df.Ks %>% filter(Helix == "F")
df.Ks.Monovalent = df.Ks %>%
  filter(Condition == "Monovalent") %>%
  filter(SE.lnK <= K_error) %>%
  filter(K <= 1/(10*10^-9)) %>%
  filter(K >= 1/(500*10^-9))
df.Ks.NTPCM = df.Ks %>%
  filter(Condition == "NTPCM") %>%
  filter(SE.lnK <= K_error) %>%
  filter(K <= 1/(10*10^-9)) %>%
  filter(K >= 1/(500*10^-9))
df.Ks.WMCM = df.Ks %>%
  filter(Condition == "WMCM") %>%
  filter(SE.lnK <= K_error) %>%
  filter(K <= 1/(20*10^-9)) %>%
  filter(K >= 1/(750*10^-9))
df.Ks.Ecoli80 = df.Ks %>%
  filter(Condition == "Ecoli80") %>%
  filter(SE.lnK <= K_error) %>%
  filter(K <= 1/(20*10^-9)) %>%
  filter(K >= 1/(750*10^-9))

vh.plot = ggplot(bind_rows(df.Ks.Monovalent, df.Ks.NTPCM, df.Ks.WMCM, df.Ks.Ecoli80),
       aes(x = invT, y = lnK,
                  ymin = lnK - SE.lnK,
                  ymax = lnK +SE.lnK,
                  color = Condition,
           shape = Condition)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  geom_pointrange() +
  ylab("ln[K (M)]") +
  xlab("1/Temperature (K)") +
  theme_classic() +
  xlim(0.003, 0.003125) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = c(0.9, 0.3))


vh.plot

####Make thermo plot####

list.files("Figures/Figure_3")

df.vh = read.csv("Figures/Figure_3/Fits_summary_vh.csv")

df.vh$Condition = factor(df.vh$Condition,
                          levels = c("Monovalent", "NTPCM", "WMCM", "Ecoli80"))


dG.plot = ggplot(data = df.vh %>% filter(Method == "1 VH plot"),
       mapping = aes(x = Condition, y = G, ymin = G - SE.G, ymax = G + SE.G)) +
  facet_wrap(~Helix, nrow = 1) +
  geom_bar(stat="identity") +
  geom_errorbar() +
  theme_classic() +
  ylab("dG 37\u00b0C (kcal/mol)") +
  scale_y_continuous(trans = 'reverse') +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 14))

dG.plot

####Read in figure 3A####

library(svglite)
library(magick)
library(rsvg)
library(grobblR)
library(grid)

list.files("Figures/Figure_3")

bitmap <- rsvg_raw('Figures/Figure_3/Figure_3A.svg', width = 600)
str(bitmap)
image <- ggdraw() + draw_image(bitmap)


####Make figure 3####

Figure_3ABC = plot_grid(image, raw.plot, vh.plot, nrow = 1,
                        labels = c("A", "B", "C"))


Figure_ABCDE = plot_grid(Figure_3ABC,
                         Tm.plot,
                         dG.plot,
                         ncol = 1,
                         labels = c("", "D", "E"))

Figure_ABCDE

ggsave("Figures/Figure_3/Figure_3.png",
       Figure_ABCDE,
       scale = 2,
       width = 7, height = 6, units = "in")
