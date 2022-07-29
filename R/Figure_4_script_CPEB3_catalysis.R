library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(viridis)

####Read in data####

list.files("Figures/Figure_4_CPEB3_catalysis")

df = read.csv("Figures/Figure_4_CPEB3_catalysis/kineticsdatafinal.csv")

head(df)

####Calculate the fraction that was cleaved####

df$f.cleaved = (df$Counts.clvd/13)/((df$Counts.clvd/21) + (df$Counts.pre/13))

hist(df$f.cleaved)

####Plot data####


df$Replicate = factor(df$Replicate)
levels(df$Replicate)

df$Condition = factor(df$Condition)
levels(df$Condition)

df$Condition = factor(df$Condition,
                      levels = c("2 mM F ",   "NTPCM","WMCM ", "Eco 80",  "25 mM F "),
                      labels = c("2 mM Free", "NTPCM",
                                 "WMCM", "Eco80", "25 mM Free"))

head(df)

ggplot(df, aes(x = Time.min, y = f.cleaved, color = Replicate)) +
  geom_point() +
  facet_wrap(~Condition, nrow = 1) +
  theme_classic() +
  xlab("Time (min)") +
  ylab("Fraction cleaved")


####Filter out long 25 mM Mg Time points####

df.25  = df %>% filter(Condition == "25 mM Free") %>% filter(Time.min <= 500)
df.no.25 = df %>% filter(Condition != "25 mM Free")

df = bind_rows(df.25, df.no.25)

####Fit by groups####

list.df.result = {}
list.df.model = {}

for (i in 1:length(unique(df$Replicate))){
  list.df.result.con = {}
  list.df.model.con = {}
  for (j in 1:length(unique(df$Condition))){
    print(unique(df$Replicate)[i])
    print(unique(df$Condition)[j])
    df.rep = df %>%
      filter(Replicate == unique(df$Replicate)[i]) %>%
      filter(Condition == unique(df$Condition)[j])
    fit = nls(f.cleaved ~ A + B*exp(-k*Time.min),
                         df.rep,
                         list(A = 0.5,
                              B = -0.5,
                              k = 0.009),
                         algorithm = "port",
                         trace = TRUE,
                         lower = c(0.1, -1, 0.0000001),
                         upper = c(1, 0, 1),
                         control = nls.control(warnOnly = TRUE))
    df.model = data.frame("Replicate" = df.rep$Replicate[1],
                          "Condition" = df.rep$Condition[1],
                          "TP" = NA,
                          "Counts.pre" = NA,
                          "Counts.clvd"= NA,
                          "Time.min" = seq(range(df.rep$Time.min)[1], range(df.rep$Time.min)[2], length.out = 100),
                          f.cleaved = NA)
    df.model$f.cleaved = predict(fit, df.model)
    list.df.model.con[[j]] = df.model
    df.result = data.frame(coef(summary(fit)))
    df.result$Parameter = c("A", "B", "k")
    df.result$Replicate = df.rep$Replicate[1]
    df.result$Condition = df.rep$Condition[1]
    list.df.result.con[[j]] = df.result
  }

  list.df.result[[i]] = bind_rows(list.df.result.con)
  list.df.model[[i]] = bind_rows(list.df.model.con)

}

df.result = bind_rows(list.df.result)
df.model = bind_rows(list.df.model)

head(df.result)

df.result$Condition = factor(df.result$Condition)

levels(df.result$Condition)

df.result$Condition = factor(df.result$Condition,
                             levels = c("2 mM Free",  "Eco80" ,     "NTPCM",      "WMCM",  "25 mM Free"))

df$f.cleaved = (df$Counts.clvd/13)/((df$Counts.clvd/21) + (df$Counts.pre/13))

hist(df$f.cleaved)

####Plot data####

df = read.csv("Figures/Figure_4_CPEB3_catalysis/kineticsdatafinal.csv")

df$f.cleaved = (df$Counts.clvd/13)/((df$Counts.clvd/21) + (df$Counts.pre/13))

df$Replicate = factor(df$Replicate)
df$Replicate = factor(df$Replicate,
                      levels = c("J1", "J2", "L1", "L2"),
                      labels = c("1", "2", "3", "4"))

df$Condition = factor(df$Condition)
levels(df$Condition)

df$Condition = factor(df$Condition,
                      levels = c("2 mM F ", "Eco 80",  "NTPCM","WMCM ",   "25 mM F "),
                      labels = c("2 mM Free", "Eco80", "NTPCM",
                                 "WMCM",  "25 mM Free"))

####Calculate k.rel####

df.result = df.result %>% filter(Parameter == "k")

df.krel = aggregate(df.result$Estimate, by = list(df.result$Condition), FUN = mean)

colnames(df.krel) = c("Condition", "k")

df.krel$SD = aggregate(df.result$Estimate, by = list(df.result$Condition), FUN = sd)[[2]]

df.krel

df.krel$k.rel = df.krel$k/df.krel$k[1]
df.krel$Y = 0.05

####C####

C = ggplot(df.result, aes(x = Condition, y = Estimate)) +
  geom_bar(position = "dodge",
           stat = "summary",
           fun = "mean") +
  geom_beeswarm(mapping = aes(shape = Replicate)) +
  geom_text(data = df.krel,
            mapping = aes(x = Condition, y = Y, label = round(k.rel, digits = 1)),
            size = 6) +
  ylab("Rate constant k (1/min)") +
  theme_classic() +
  ylim(0, 0.05) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        legend.position = "none")

C

####B####

B = ggplot(mapping =  aes(x = Time.min,
                           shape = Replicate,
                           group = Replicate,
                           y = f.cleaved,
                           color = Condition)) +
  geom_point(data = df) +
  geom_line(data = df.model) +
  facet_wrap(~Condition, ncol = 1) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  scale_x_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(size = 16),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  xlab("Time (min)") +
  ylab("Fraction cleaved")

B

ggplot(mapping =  aes(x = Time.min,
                      shape = Replicate,
                      group = Replicate,
                      y = f.cleaved,
                      color = Condition)) +
  geom_point(data = df) +
  geom_line(data = df.model) +
  facet_wrap(~Condition, ncol = 1) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme_classic() +
  #scale_x_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(size = 16),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  xlab("Time (min)") +
  ylab("Fraction cleaved")

ggsave("Figures/Figure_4_CPEB3_catalysis/Figure_4B_CPEB3_catalysis_no_log.png",
       width = 1.5, height = 4, scale = 2.5, bg = "white")



####A###

library(svglite)
library(magick)
library(rsvg)
library(grobblR)
library(grid)

list.files("Figures/Figure_4_CPEB3_catalysis/")

bitmap <- rsvg_raw('Figures/Figure_4_CPEB3_catalysis/CPEB3_secondary_structure.svg', width = 600)
A <- ggdraw() + draw_image(bitmap)

A

####Figure 4D####

####Load in data####

list.files("Figures/Figure_4_CPEB3_catalysis")

E.coli <- read.csv("Figures/Figure_1/Top_15_E.coli_metabolites_edited.csv")

head(E.coli)

Yeast = read.csv("Figures/Figure_4_CPEB3_catalysis/Table_1_final.csv") %>% filter(Species == "Yeast")
Mammal = read.csv("Figures/Figure_4_CPEB3_catalysis/Table_1_final.csv") %>% filter(Species == "Mammal iBMK")

####Sum metabolites####

#E.coli

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

#Yeast

Yeast$Mg.binding.strength = NA

Yeast$Mg.binding.strength[which(Yeast$Kd.app >= 2)] = "weak"
Yeast$Mg.binding.strength[which(Yeast$Kd.app <= 2)] = "strong"

Sum.concentration = c()

for (i in 1:length(Yeast$Metabolites)){
  Sum.concentration[i] = sum(Yeast$Concentration[1:i])
}

Yeast$Sum.concentration = Sum.concentration

Sum.concentration[80]*0.80

Yeast = Yeast %>% filter(Sum.concentration <= 193) %>% select(Species, Metabolites, Concentration, Mg.binding.strength)

sum(Yeast$Concentration)

table(Yeast$Mg.binding.strength)

Yeast.weak = Yeast %>% filter(Mg.binding.strength == "weak")

Sum.concentration.weak = c()

for (i in 1:length(Yeast.weak$Metabolites)){
  Sum.concentration.weak[i] = sum(Yeast.weak$Concentration[1:i])
}

Yeast.strong = Yeast %>% filter(Mg.binding.strength == "strong")

Sum.concentration.strong = c()

for (i in 1:length(Yeast.strong$Metabolites)){
  Sum.concentration.strong[i] = sum(sum(Yeast.weak$Concentration), Yeast.strong$Concentration[1:i])
}

Sum.concentration = c(Sum.concentration.weak)

Yeast = bind_rows(Yeast.weak, Yeast.strong)

Yeast$Sum.concentration = Sum.concentration

#Mammal

Mammal$Mg.binding.strength = NA

Mammal$Mg.binding.strength[which(Mammal$Kd.app >= 2)] = "weak"
Mammal$Mg.binding.strength[which(Mammal$Kd.app <= 2)] = "strong"

Sum.concentration = c()

for (i in 1:length(Mammal$Metabolites)){
  Sum.concentration[i] = sum(Mammal$Concentration[1:i])
}

Mammal$Sum.concentration = Sum.concentration

Sum.concentration[90]*0.80

Mammal = Mammal %>% filter(Sum.concentration <= 141) %>% select(Species, Metabolites, Concentration, Mg.binding.strength)

sum(Mammal$Concentration)

table(Mammal$Mg.binding.strength)

Mammal.weak = Mammal %>% filter(Mg.binding.strength == "weak")

Sum.concentration.weak = c()

for (i in 1:length(Mammal.weak$Metabolites)){
  Sum.concentration.weak[i] = sum(Mammal.weak$Concentration[1:i])
}

Mammal.strong = Mammal %>% filter(Mg.binding.strength == "strong")

Sum.concentration.strong = c()

for (i in 1:length(Mammal.strong$Metabolites)){
  Sum.concentration.strong[i] = sum(sum(Mammal.weak$Concentration), Mammal.strong$Concentration[1:i])
}

Sum.concentration = c(Sum.concentration.weak, Sum.concentration.strong)

Mammal = bind_rows(Mammal.weak, Mammal.strong)

Mammal$Sum.concentration = Sum.concentration


E.coli$Species = "Eco80"
Yeast$Species = "Yeast80"
Mammal$Species = "Mammal80"

E.coli = E.coli %>% select(Species, Metabolites, Concentration, Mg.binding.strength, Sum.concentration)
head(E.coli)
head(Yeast)

df.met = bind_rows(E.coli, Yeast, Mammal)

####Make Figure 4D####

df.met$Mg.binding.strength = factor(df.met$Mg.binding.strength,
                                    levels = c("other", "strong", "weak"),
                                    labels = c("other", "NTPCM", "WMCM"))

library(viridis)

Figure_4D = ggplot(df.met, aes(x = Species, y = Concentration, fill = Mg.binding.strength, label = Metabolites)) +
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

Figure_4D

####Consolidate plots####

AC = plot_grid(A, C, ncol = 1,
               labels = c("A", "C"),
               label_size = 20,
               rel_heights = c(1, 1.3))
ABC = plot_grid(AC, B,
                labels = c("", "B"),
                label_size = 20)

ABC

ABCD = plot_grid(ABC, Figure_4D,
                 labels = c("", "D"),
                 ncol = 1,
                 rel_heights = c(3, 1),
                 label_size = 20)



#ggsave("Figures/Figure_4_CPEB3_catalysis/Figure_4_CPEB3_catalysis.svg",
#       width = 3.3, height = 4, scale = 2.5, bg = "white")

