setwd("~/Jacob/Research/Manuscripts/JPSiegMetaboMetaloRNA/Figures/SI_Figure_1")

library(tidyverse)
library(viridis)
library(cowplot)
devtools::load_all("~/Jacob/R_packages/MetaboMgITC2")

####Load in predicted binding constants####

list.files("../Figure_1")

df = read.csv("../Figure_1/SI_Table_X_ITC_fit_results.csv")

head(df)

K.app = c()
dH = c()

for (i in 1:length(df$Metabolite)){
  Kd.app = Kd.app.calc(df$Metabolite[i],
                       Temperature = df$Temp.[i],
                       constants.path = "~/Jacob/R_packages/MetaboMgITC2/Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv")
  dH[i] = Kd.app.calc(df$Metabolite[i],
                   Temperature = df$Temp.[i],
                   output.dH = TRUE,
                   constants.path = "~/Jacob/R_packages/MetaboMgITC2/Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv")
  K.app[i] = 1000/(Kd.app)
}

df$K.app = K.app
df$dH.pred = dH


head(df)

#####SI Figure 1A#####

error = 0.5

df.poly = data.frame(x = c(1,
                           1*(1+error),
                           10000,
                           10000,
                           10000*(1-error),
                           1),
                     y = c(1,
                           1,
                           10000*(1-error),
                           10000,
                           10000,
                           1*(1+error)))

SI_Figure_1A = ggplot() +
  geom_polygon(data = df.poly, mapping = aes(x = x, y = y), fill = "grey") +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = df, mapping = aes(x = K.app, y = K, shape = Metabolite, color = Temp.)) +
  scale_x_continuous(trans = "log10", limits = c(1,10000), expand = c(0,0)) +
  scale_y_continuous(trans = "log10", limits = c(1,10000), expand = c(0,0)) +
  scale_color_viridis() +
  coord_fixed() +
  theme_classic() +
  xlab("Calculated K (1/M)")+
  ylab("ITC K (1/M)")

#####SI Figure 1A#####

error = 0.5

df.poly = data.frame(x = c(1,
                           1*(1+error),
                           5,
                           5,
                           5*(1-error),
                           1),
                     y = c(1,
                           1,
                           5*(1-error),
                           5,
                           5,
                           1*(1+error)))

SI_Figure_1B = ggplot() +
  geom_polygon(data = df.poly, mapping = aes(x = x, y = y), fill = "grey") +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = df, mapping = aes(x = dH.pred, y = dH, shape = Metabolite, color = Temp.)) +
  scale_x_continuous(trans = "log10", limits = c(1,5), expand = c(0,0)) +
  scale_y_continuous(trans = "log10", limits = c(1,5), expand = c(0,0)) +
  scale_color_viridis() +
  coord_fixed() +
  theme_classic() +
  xlab("Calculated dH (kcal/mol)")+
  ylab("ITC dH (kcal/mol)")


SI_Figure_1 = plot_grid(SI_Figure_1A, SI_Figure_1B, labels = c("A", "B"))

?ggsave
ggsave("SI_Figure_1.svg", SI_Figure_1, scale = 2, width = 8, height = 3.3, units = "in")
