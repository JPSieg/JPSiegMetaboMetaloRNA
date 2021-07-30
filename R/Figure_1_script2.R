setwd("~/Jacob/Research/Manuscripts/JPSiegMetaboMetaloRNA/Figures/Figure_1")

####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(scales)
devtools::load_all("~/Jacob/R_packages/MetaboMgITC2")

options(scipen=10000)

####Load in data####

list.files()

E.coli <- read.csv("Top_15_E.coli_metabolites_edited.csv")

head(E.coli)

####Sum metabolites####

E.coli.weak = E.coli %>% filter(Mg.binding.strength == "weak")

Sum.concentration.weak = c()

for (i in 1:length(E.coli.weak$Metabolites)){
  Sum.concentration.weak[i] = sum(E.coli.weak$Concentration[1:i])
}

E.coli.strong = E.coli %>% filter(Mg.binding.strength == "strong")

Sum.concentration.strong = c()

for (i in 1:length(E.coli.strong$Metabolites)){
  Sum.concentration.strong[i] = sum(sum(E.coli.weak$Concentration), E.coli.strong$Concentration[1:i])
}

Sum.concentration = c(Sum.concentration.weak, Sum.concentration.strong)

E.coli = bind_rows(E.coli.weak, E.coli.strong)

E.coli$Sum.concentration = Sum.concentration

sum(E.coli$Concentration)

df.total =  data.frame("Metabolites" = c("103 metabolites"),
                       "Concentration" = c(243 - sum(E.coli$Concentration)),
                       "Metabolites.sum" = c(243),
                       "Kd" = c(0),
                       "Mg.binding.strength" = c("Other"),
                       "Edited" = c(NA),
                       "Sum.concentration" = c(243))

E.coli = bind_rows(E.coli, df.total)

####Make Figure 1A####

Figure_1A = ggplot(E.coli, aes(x = "", y = Concentration, fill = Mg.binding.strength, label = Metabolites)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  geom_text(mapping = aes(y = Sum.concentration), nudge_y = -3, color = "white", size = 3) +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  theme_minimal()+
  theme(axis.line.y = element_line(colour = 'black', size = 1.5),
        axis.line.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 8,
                                   angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.5, 0.2),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16))

Figure_1A

####Load in and fit ITC data####

df = read.csv("ITC_data_index.csv")

list.fits = {}
list.df.fit = {}

colnames(df)

for (i in 1:length(df$Metabolite)){
  print(paste(df$Metabolite[i], " at ", df$Temperature[i], "C", sep = ""))
  df.Cell = read.itc(c(paste("ITC_data_files", df$Cell[i], sep = "/"),
                       df$Syrringe.X[i], df$Syringe.M[i], df$Syringe.C[i],
                       df$Cell.X[i], df$Cell.M[i], df$Cell.C[i]))
  df.blank = read.itc(c(paste("ITC_data_files", df$Blank[i], sep = "/"),
                        df$Syrringe.X[i], df$Syringe.M[i], df$Syringe.C[i],
                        df$Cell.X[i], df$Cell.M[i], df$Cell.C[i]))
  list.fits[[i]] = MetaboMgITC(df.Cell, df.blank,
                               Fit.start = list(H = df$H.start[i], K = df$K.start[i]),
                               Saturation.threshold = df$Sat.threshold[i])
  list.df.fit[[i]] = data.frame(t(c(list.fits[[i]]$Table[,2], list.fits[[i]]$Table[,3])))
  colnames(list.df.fit[[i]]) = c(as.character(t(list.fits[[i]]$Table[,1])), paste("Std.error.", as.character(t(list.fits[[i]]$Table[,1])), sep = ""))
  list.df.fit[[i]]$Metabolite = df$Metabolite[i]
}

df.final = bind_rows(list.df.fit)

df.final$Kd = 1000/df.final$K

Std.error.Kd = c()

for (i in 1:length(df.final$n)){
  K = df.final$K[i]
  dK = df.final$Std.error.K[i]
  Kd = z ~ 1000/K
  Std.error.Kd[i] = abs(dK*eval(D(Kd[[3]], "K")))
}


df.final$Std.error.Kd = Std.error.Kd

colnames(df.final)

df.final = df.final %>% select(Metabolite, Temp., Std.error.Temp.,
                               K, Std.error.K,
                               Kd, Std.error.Kd,
                               n, Std.error.n,
                               dG, Std.error.dG,
                               dH, Std.error.dH,
                               dS, Std.error.dS,
                               Saturation, Std.error.Saturation)
df.final

####Write SI table####

write.csv(df.final, "SI_Table_X_ITC_fit_results.csv", row.names = FALSE)

Figure_1B = list.fits[[9]]$Plot

####Make Figure 1C####

df.final$ymin = df.final$K - df.final$Std.error.K
df.final$ymax = df.final$K + df.final$Std.error.K
df.final$xmin = df.final$Temp. - df.final$Std.error.Temp.
df.final$xmax = df.final$Temp. + df.final$Std.error.Temp.


df.rect = data.frame(x = c(24, 51, 51, 24), y = c(10, 10, 1/0.002, 1/0.002))

head(df.final)

one_over_trans = function() trans_new("one_over", function(x) 1/(x + 273.15), function(x)((1/x)  - 273.15))

Figure_1C = ggplot() +
  geom_polygon(data = df.rect, mapping = aes(x = x, y = y), fill = "grey") +
  geom_hline(yintercept = 1/0.002, size = 1.5) +
  stat_smooth(data = df.final, mapping = aes(x = Temp.,
                                             y = K,
                                             color = Metabolite,), method = "lm", se = FALSE) +
  geom_pointrange(data = df.final, mapping = aes(x = Temp., color = Metabolite, shape = Metabolite, y = K, xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), size = 1) +
  annotate("text", y = c(650, 403, 12), x = c(45, 45, 30), label = c("Strong", "Weak", "transformed = 1/(x+273.15)"), size = 5) +
  scale_color_manual(values = viridis(9)) +
  scale_y_continuous(trans = "log10", limits = c(10, 15000), expand = c(0, 0)) +
  scale_x_continuous(trans = "one_over", expand = c(0, 0)) +
  ylab("K (1/M)") +
  xlab("Temperature (\u00b0C)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 12),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.55))

Figure_1C

####Predict how strong amino acids bind to Mg at different temperatures####

#?Kd.app.calc

list.files()

df = read.csv("Top_15_E.coli_amino_acids.csv")

head(df)

Kd.app = c()

for (i in 1:length(df$Metabolites)){
  Kd.app[i] = Kd.app.calc(df$Metabolites[i],
              Temperature = df$Temperature[i],
              constants.path = "~/Jacob/R_packages/MetaboMgITC2/Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv")
}

df$K.app = 1000/Kd.app

####Make Figure 1D####

length(unique(df$Metabolites))

df.rect = data.frame(x = c(24, 51, 51, 24), y = c(0.001, 0.001, 150, 150))

Figure_1D = ggplot() +
  geom_polygon(data = df.rect, mapping = aes(x = x, y = y), fill = "grey") +
  geom_line(data = df, mapping = aes(x = Temperature, color =  Metabolites, y = K.app), size = 2) +
  scale_color_manual(values = viridis(7)) +
  scale_y_continuous(trans = "log10", limits = c(0.001, 150), expand = c(0, 0)) +
  scale_x_continuous(trans = "one_over", expand = c(0, 0)) +
  annotate("text", y = c(650, 403, 0.0015), x = c(45, 45, 30), label = c("Strong", "Weak", "transformed = 1/(x+273.15)"), size = 5) +
  ylab("K (1/M)") +
  xlab("Temperature (\u00b0C)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 12),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.5))

Figure_1D

write.csv(df, "SI_Table_X_predicted_binding_constants.csv", row.names = FALSE)

####Put together final figure####

#?plot_grid

Figure_1CD = plot_grid(Figure_1C, Figure_1D, nrow = 2, labels = c("C", "D"), label_size = 20)
Figure_1BCD = plot_grid(Figure_1B, Figure_1CD, rel_widths = c(1, 1.2), labels = c("B"), label_size = 20)
Figure_1 = plot_grid(Figure_1A, Figure_1BCD, rel_widths = c(1, 5), labels = c("A"), label_size = 20)

ggsave("Figure_1.svg", Figure_1, scale = 5, units = "in", width = 3.3, height = 2.7, bg = "white")
