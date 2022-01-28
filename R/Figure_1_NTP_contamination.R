library(tidyverse)
library(ggbeeswarm)
library(viridis)

df.HQS = read.csv("Figures/Figure_1/HQS_data.csv")

df = df.HQS %>% filter(Metabolites == "NTPCM")

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10))

df$I.norm = (df$Emission - coef(fit)[2])/(coef(fit)[1]- coef(fit)[2])

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Figure_free.Mg_EDTA_check = ggplot(df %>% filter(Conc.Mg %in% c(0, 1.5228)), aes(x = Sample, y = Mg.free, color = EDTA)) +
  geom_beeswarm() +
  scale_color_manual(values = viridis(5)) +
  ylab("[Mg] free (mM)") +
  xlab("") +
  theme_classic() +
  ggtitle("No Mg2+ control") +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = c(0.75, 0.75))

E.coli <- read.csv("Figures/Figure_1/Top_15_E.coli_metabolites_edited.csv") %>%
  filter(Mg.binding.strength == "strong")

df.contam = df %>% filter(Conc.Mg == 0) %>%
  filter(EDTA == "EDTA = 0 mM") %>%
  filter(Sample == "No chelator")

Contam = c()

for (i in 1:length(df.contam$Mg.free)){
  Contam[i] = df.contam$Mg.free[i] + sum((E.coli$Concentration*df.contam$Mg.free[i])/(df.contam$Mg.free[i] + E.coli$Kd))
}

print(mean(Contam))
