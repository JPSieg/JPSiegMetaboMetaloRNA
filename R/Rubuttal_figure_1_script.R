library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(viridis)

####Load in data####

list.files("Figures/Figure_4_CPEB3_catalysis")

df.Eco80 <- read.csv("Figures/Figure_1/Top_15_E.coli_metabolites_edited.csv")

df.Eco80

df = read.csv("Figures/Figure_4_CPEB3_catalysis/Table_1_final.csv") %>%
  filter(Species == "E.coli") %>%
  filter(!Metabolites %in% df.Eco80$Metabolites)

####Sum metabolites####

#E.coli

df$Mg.binding.strength = NA
df$Mg.binding.strength[which(df$Kd.app <= 2)] = "strong"
df$Mg.binding.strength[-which(df$Kd.app <= 2)] = "weak"

df.weak = df %>% filter(Mg.binding.strength == "weak")

Sum.concentration.weak = c()

for (i in 1:length(df.weak$Metabolites)){
  Sum.concentration.weak[i] = sum(df.weak$Concentration[1:i])
}

df.strong = df %>% filter(Mg.binding.strength == "strong")

Sum.concentration.strong = c()

for (i in 1:length(df.strong$Metabolites)){
  Sum.concentration.strong[i] = sum(sum(df.strong$Concentration), df.strong$Concentration[1:i])
}

Sum.concentration = c(Sum.concentration.weak, Sum.concentration.strong)

df = bind_rows(df.weak, df.strong)

df$Sum.concentration = Sum.concentration

sum(df$Concentration)

Total = 194.88 + sum(df$Concentration)

df$Mg.binding.strength = factor(df$Mg.binding.strength,
                                    levels = c("strong", "weak"),
                                    labels = c("NTPCM", "WMCM"))

library(viridis)

ggplot(df, aes(x = "", y = Concentration, fill = Mg.binding.strength, label = Metabolites)) +
  geom_bar(width = 0.8, stat = "identity", color = "black", size = 1) +
  scale_fill_manual(values = viridis(n =  7)[c(7, 3, 1)]) +
  theme_classic()+
  #scale_x_discrete(limits = factor("",)) +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 10),
        legend.background = element_blank()) +
  ylab("[Metabolites] (mM)") +
  coord_flip()

list.files("Reviewer_response_figures")

#ggsave("Reviewer_response_figures/R_figure_1.svg", width = 5.3, height = 1.66)

write.csv(df, "Reviewer_response_figures/Other_metabolites.csv", row.names = FALSE)
sum(df.weak$Concentration)/sum(df$Concentration)
sum(df.strong$Concentration)/sum(df$Concentration)

sum(df.Eco80$Concentration[which(df.Eco80$Mg.binding.strength == "strong")])/194.88
sum(df.Eco80$Concentration[which(df.Eco80$Mg.binding.strength == "weak")])/194.88

df$Concentration/Total

