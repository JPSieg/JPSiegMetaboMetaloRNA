
####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(scales)
library(MuMIn)
#devtools::load_all("~/Jacob/R_packages/MetaboMgITC2")

options(scipen=10000)

####Functions####
find.Mg.total = function(fit,
                         Mg.total.start = c(1:100),
                         Free.Mg.mM = 2.75){
  terms = c()

  for(i in length(fit$coefficients):1){
    terms[i] = paste(fit$coefficients[i], "*(Mg.total)^(", i - 1, ")", sep = "")
    if (i == length(fit$coefficients)){
      formula.text = paste(terms[i], sep = "")
    }else{
      formula.text = paste(formula.text, terms[i], sep = "+")
    }
  }

  first.recursion = TRUE
  Go = TRUE
  while(Go){
    if (first.recursion){
      Mg.free = c()
      for (i in 1:length(Mg.total.start)){
        Mg.total = Mg.total.start[i]
        formula.text.rec = gsub("Mg.total", Mg.total, formula.text)
        formula.text.rec = gsub("NA", 0, formula.text.rec)
        Mg.free[i] = eval(parse(text = formula.text.rec))
        first.recursion = FALSE
      }
    }else{
      Mg.free.error = abs(Mg.free - Free.Mg.mM)
      if(min(Mg.free.error) <= 0.0001){
        output = Mg.total.start[which.min(Mg.free.error)]
        Go = FALSE
      }else{
        best.Mg.total = Mg.total.start[which.min(Mg.free.error)]
        Mg.total.start = seq(best.Mg.total - min(Mg.free.error), best.Mg.total + min(Mg.free.error), length.out = 100)
        Mg.free = c()
        for (i in 1:length(Mg.total.start)){
          Mg.total = Mg.total.start[i]
          formula.text.rec = gsub("Mg.total", Mg.total, formula.text)
          formula.text.rec = gsub("NA", 0, formula.text.rec)
          Mg.free[i] = eval(parse(text = formula.text.rec))
          first.recursion = FALSE
        }
      }
    }
  }
  print(output)
  output = output
}



####Load in data####

list.files("Figures")

E.coli <- read.csv("Figures/Figure_1/Top_15_E.coli_metabolites_edited.csv")

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
                       "Mg.binding.strength" = c("other"),
                       "Edited" = c(NA),
                       "Sum.concentration" = c(243))

E.coli = bind_rows(E.coli, df.total)

####Make Figure 1A####

E.coli$Mg.binding.strength = factor(E.coli$Mg.binding.strength,
       levels = c("other", "strong", "weak"),
       labels = c("other", "NTPCM", "WMCM"))

Figure_1A = ggplot(E.coli, aes(x = "", y = Concentration, fill = Mg.binding.strength, label = Metabolites)) +
  geom_bar(width = 0.8, stat = "identity", color = "black") +
  annotate("line", x = c(1.4,1.4), y = c(0,195), size = 2, color = viridis(n = 7)[6]) +
  annotate("line", x = c(1.4,1.4), y = c(195,243), size = 2, color = "black") +
  annotate("text", x = 1.53, y = 100, label = "15 metabolites = 80% E. coli metabolome = Ecoli80", size = 4,
           color = viridis(n = 7)[6]) +
  annotate("text", x = 1, y = 220, label = "10%\n228 other\nmetabolites", size = 4,
           color = "black") +
  scale_fill_manual(values = viridis(n =  7)[c(7, 3, 1)]) +
  theme_classic()+
  #scale_x_discrete(limits = factor("",)) +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(color = "Black", size = 8,
                                   angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.15, 0.55),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 10),
        legend.background = element_blank()) +
  ylab("[Metabolites] (mM)") +
  coord_flip()

Figure_1A

#####Figure 1 B-G####


list.files("Figures/Figure_1")

analyze.HQS = function(df = df.HQS %>% filter(Metabolites == "WMCM"),
                       df.model = df.model.WMCM,
                       color = viridis(n =  7)[1],
                       Labels = c("B", "E")){

  df = df %>% filter(EDTA == "EDTA = 0 mM")

  fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
            df %>% filter(Sample == "No chelator"),
            start = list(I.max = 150000, I.min = 0, K = 10))

  fit.form = function(Conc.Mg){
    Emission = (coef(fit)[1] - coef(fit)[2])*(coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg)) + coef(fit)[2]
  }

  df$I.norm = (df$Emission - coef(fit)[2])/(coef(fit)[1]- coef(fit)[2])

  fit1 = fit
  fit.form = function(Conc.Mg){
    Emission = (coef(fit1)[3]*Conc.Mg/(1 + coef(fit1)[3]*Conc.Mg))
  }

  Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color = Sample)) +
    geom_point() +
    theme_classic() +
    geom_function(fun = fit.form, color = "dimgrey") +
    ylab("HQS emmission") +
    xlab("[Mg] total (mM)") +
    ggtitle(df$Metabolites[1]) +
    scale_color_manual(values = c(color, "dimgrey")) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.text = element_text(color = "Black", size = 10),
          legend.title = element_text(color = "Black", size = 16),
          legend.position = c(0.8, 0.3),
          plot.title = element_text(color = "Black", size = 14,hjust = 0.5))

  df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

  list.fit = {}

  for (i in 1:10){
    list.fit[[i]] = lm( as.vector(Mg.free) ~  poly(Conc.Mg, i, raw=TRUE), df %>% filter(Sample == "Chelator"))
  }

  df.model.sel = data.frame(model.sel(list.fit))

  best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
  best.polynomial = list.fit[[best.polynomial.order]]

  df.no_edta = df %>% filter(Sample == "Chelator")  %>% filter(EDTA == "EDTA = 0 mM")

  df.no_edta$model = predict(list.fit[[best.polynomial.order]])

  df.no_edta2 = df %>% filter(Sample == "No chelator")  %>% filter(EDTA == "EDTA = 0 mM")
  df.no_edta2$model = NA

  df.no_edta = bind_rows(df.no_edta, df.no_edta2)

  free.Mg = 2

  Mg.total = find.Mg.total(best.polynomial, Free.Mg.mM = free.Mg)

  fit = best.polynomial

  terms = c()

  for(i in length(fit$coefficients):1){
    terms[i] = paste(fit$coefficients[i], "*(Total.Mg)^(", i - 1, ")", sep = "")
    if (i == length(fit$coefficients)){
      formula.text = paste(terms[i], sep = "")
    }else{
      formula.text = paste(formula.text, terms[i], sep = "+")
    }
  }

  Figure_Mg_free = ggplot() +
    geom_hex(data = df.model, mapping = aes(x = Conc.Mg, y = Mg.free), bins = 100) +
    geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
    geom_line(data = df.no_edta %>% filter(!is.na(model)), mapping = aes(x = Conc.Mg, y = model), color = color) +
    geom_point(data = df.no_edta, mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
    theme_classic() +
    scale_fill_viridis(option = "rocket") +
    geom_hline(yintercept = 2,
                 color = "red") +
    #annotate("text", x = 30, y = 0.01, label = paste(round(Mg.total, digits = 2), " mM total Mg2+"), color = "red") +
    annotate("text", x = 0.7, y = 10, label = paste(round(free.Mg, digits = 2), " mM free Mg2+"), color = "red") +
    scale_color_manual(values = c(color, "dimgrey")) +
    scale_y_continuous(trans = "log10",
                       breaks = c(0.1, 1, 10, 100, 200),
                       labels = c("0.1", "1.0", "10", "100", "200"),
                       lim = c(0.1, 250)) +
    scale_x_continuous(trans = "log10",
                       breaks = c(0.1, 1, 10, 100, 200),
                       labels = c("0.1", "1.0", "10", "100", "200"),
                       lim = c(0.1, 250)) +
    ylab("[Mg] free (mM)") +
    xlab("[Mg] total (mM)") +
    ggtitle(df$Metabolites[1]) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.position = "none",
          plot.title = element_text(color = "Black", size = 14,hjust = 0.5))

  output = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, labels = Labels)

}

df.HQS = read.csv("Figures/Figure_1/HQS_data.csv")

df.AC.model = read.csv("Figures/Figure_1/Modeled_AC_MCM_concentrations.csv")

#NTPCM
df.model.NTPCM = df.AC.model %>% select(Mg.T, Mg.free.NTP)
colnames(df.model.NTPCM) = c('Conc.Mg', "Mg.free")

Figure_1BE = analyze.HQS(df.HQS %>% filter(Metabolites == "NTPCM"),
                         df.model.NTPCM,
                         viridis(n =  7)[3],
                         Labels = c("B", "E"))

#WMCM
df.model.WMCM = df.AC.model %>% select(Mg.T, Mg.free.WMCM)
colnames(df.model.WMCM) = c('Conc.Mg', "Mg.free")

Figure_1CF = analyze.HQS(df.HQS %>% filter(Metabolites == "WMCM"),
                         df.model.WMCM,
                         viridis(n =  7)[1],
                         Labels = c("C", "F"))

#Ecoli80
df.model.Ecoli80 = df.AC.model %>% select(Mg.T, Mg.free)
colnames(df.model.Ecoli80) = c('Conc.Mg', "Mg.free")

Figure_1DG = analyze.HQS(df.HQS %>% filter(Metabolites == "Ecoli80"),
                         df.model.Ecoli80,
                         viridis(n =  7)[6],
                         Labels = c("D", "G"))


####Figure 1 H####

df.conc.m = read.csv("Figures/Figure_1/Modeled_AC_metabolite_concentrations.csv")

df.conc.m$MCM = factor(df.conc.m$MCM,
                   levels = c("NTPCM","WMCM", "Ecoli80"),
                   labels = c("Strong", "Weak", "Ecoli80"))


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

Figure_1H = ggplot(df.conc.m, aes(x = MCM,
                      y = MCM.conc,
                      fill = MCM)) +
  geom_violin() +
  stat_summary(fun.data = quantiles_95, geom="boxplot", alpha = 0) +
  scale_fill_manual(values = viridis(n =  7)[c(3,1,6)]) +
  theme_classic() +
  ylim(0, 250) +
  theme(axis.text = element_text(color = "Black", size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "Black", size = 16)) +
  ylab("[Metabolites] (mM)") +
  coord_flip()

####Make Figure_1I####

head(df.AC.model)

colnames(df.AC.model)[1] = "Iteration"

MCM = c()
Mg.T = c()

for (i in 1:length(unique(df.AC.model$Iteration))){
  df = df.AC.model %>% filter(Iteration == unique(df.AC.model$Iteration)[i])
  MCM = c(MCM, "NTPCM")
  Mg.T = c(Mg.T, df$Mg.T[which.min(abs(2 - df$Mg.free.NTP))])
  MCM = c(MCM, "WMCM")
  Mg.T = c(Mg.T, df$Mg.T[which.min(abs(2 - df$Mg.free.WMCM))])
  MCM = c(MCM, "Ecoli80")
  Mg.T = c(Mg.T, df$Mg.T[which.min(abs(2 - df$Mg.free))])
}


df.Mg.T = data.frame(MCM, Mg.T)

df.Mg.T$MCM = factor(df.Mg.T$MCM,
                       levels = c("NTPCM","WMCM", "Ecoli80"),
                       labels = c("Strong", "Weak", "Ecoli80"))

Figure_1I = ggplot(df.conc.m, aes(x = MCM,
                                  y = Mg.T,
                                  fill = MCM)) +
  geom_violin() +
  stat_summary(fun.data = quantiles_95, geom="boxplot", alpha = 0) +
  geom_hline(yintercept = c(6.44, 24.94, 31.02), color = "red") +
  scale_fill_manual(values = viridis(n =  7)[c(3,1,6)]) +
  theme_classic() +
  theme(axis.text = element_text(color = "Black", size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "Black", size = 16)) +
  ylab("[Mg Total] (mM) for 2 mM Mg Free") +
  coord_flip()


Figure_1BCDEFG = plot_grid(Figure_1BE, Figure_1CF, Figure_1DG, nrow = 1)

Figure_1HI = plot_grid(Figure_1H, Figure_1I, nrow = 1, labels = c("H", "I"))

Figure_1ABCDEFGHI = plot_grid(Figure_1A, Figure_1BCDEFG, Figure_1HI, labels = "A", ncol = 1, rel_heights = c(0.5,2,1))

ggsave("Figures/Figure_1/Figure_1ABCDEFGHI.png", Figure_1ABCDEFGHI, width = 3.3, height = 4, units = "in", scale = 3)
