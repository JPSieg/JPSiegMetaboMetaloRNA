
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
  geom_bar(width = 0.8, stat = "identity", color = "black", size = 1) +
  annotate("rect", xmin = 0.6, xmax = 1.4, ymin = 0, ymax = 195,
           size = 2, color = "black", alpha = 0.1) +
  annotate("text", x = 1.5, y = 100, label = "15 metabolites = 80% E. coli metabolome = Eco80", size = 4,
           color = "black") +
  annotate("text", x = 1, y = 220, label = "228\nmetabolites\n = other 20%", size = 4,
           color = "black") +
  scale_fill_manual(values = viridis(n =  7)[c(7, 3, 1)]) +
  theme_classic()+
  #scale_x_discrete(limits = factor("",)) +
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(color = "Black", size = 8,
                                   angle = 45, hjust = 1, vjust = 1),
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

analyze.HQS = function(df = df.HQS %>% filter(Metabolites == "Eco80"),
                       df.model = df.model.Ecoli80,
                       color = viridis(n =  7)[1],
                       Labels = c("B", "E"),
                       xlimits = c(0, 75),
                       ylimits = c(-1, 20)){


  df = df %>% filter(EDTA == "EDTA = 0 mM")

  fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
            df %>% filter(Sample == "No chelator"),
            start = list(I.max = 150000, I.min = 0, K = 10),
            trace = TRUE)

  output.HQS = data.frame(coef(summary(fit)))
  output.HQS$Condition = df$Metabolites[1]

  fit.form = function(Conc.Mg){
    Emission = (coef(fit)[1] - coef(fit)[2])*(coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg)) + coef(fit)[2]
  }

  df$I.norm = (df$Emission - coef(fit)[2])/(coef(fit)[1]- coef(fit)[2])

  fit1 = fit
  fit.form = function(Conc.Mg){
    Emission = (coef(fit1)[3]*Conc.Mg/(1 + coef(fit1)[3]*Conc.Mg))
  }

  df$Sample = factor(df$Sample,
                     levels = c("No chelator", "Chelator"))

  Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, shape = Sample, color = Sample)) +
    geom_point() +
    theme_classic() +
    geom_function(fun = fit.form, color = "dimgrey") +
    ylab("HQS emission") +
    xlab(expression("Total [Mg^2+] (mM)")) +
    ggtitle(df$Metabolites[1]) +
    scale_color_manual(values = c("dimgrey", "black")) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 210)) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.text = element_text(color = "Black", size = 10),
          legend.title = element_blank(),
          legend.position = c(0.30, 0.83),
          plot.title = element_text(color = "Black", size = 14,hjust = 0.5))

  df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

  list.fit = {}

  for (i in 1:10){
    list.fit[[i]] = lm(Mg.free ~  poly(Conc.Mg, i, raw=TRUE), df %>% filter(Sample == "Chelator"))
  }

  df.model.sel = data.frame(model.sel(list.fit))

  best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
  best.polynomial = list.fit[[best.polynomial.order]]

  df.no_edta = df %>% filter(Sample == "Chelator")  %>% filter(EDTA == "EDTA = 0 mM")

  df.no_edta$model = predict(list.fit[[best.polynomial.order]])

  coef(summary(list.fit[[best.polynomial.order]]))

  plot(df.no_edta$Conc.Mg, df.no_edta$model)

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
  formula.text = gsub("Total.Mg", "x", formula.text)

  lm.fun = function(x){
    eval(parse(text = formula.text))
  }

  colnames(df.model) = c("Conc.Mg", "Mg.free")

  Figure_Mg_free = ggplot() +
    geom_hex(data = df.model %>% filter(Conc.Mg < 200), mapping = aes(x = Conc.Mg, y = Mg.free), bins = 50) +
    geom_polygon(data = NULL, mapping = aes(x = c(0.1, 200, 200, 0.1), y = c(0.5, 0.5, 3, 3)),
                 fill = "red", alpha = 0.25) +
    annotate("segment", y = 2, yend = 2, x = 0.1, xend = 200,
             color = "red") +
    geom_point(data = df.no_edta %>% filter(Sample == "Chelator"), mapping = aes(x = Conc.Mg, y = Mg.free), shape = "triangle") +
    theme_classic() +
    geom_line(data = df.no_edta, mapping = aes(x = Conc.Mg, y = model),  color = "black") +
    scale_fill_viridis(option = "D") +
    #annotate("text", x = 30, y = 0.01, label = paste(round(Mg.total, digits = 2), " mM total Mg2+"), color = "red") +
    #annotate("text", x = 40, y = 5, label = paste(round(free.Mg, digits = 2), " mM free Mg2+"), color = "red") +
    ylab(bquote('Free'~'Mg'^'2+'~'(mM)')) +
    xlab(bquote('Total'~'Mg'^'2+'~'(mM)')) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 210)) +
    scale_y_continuous(lim = ylimits) +
    geom_vline(xintercept = c(40, 200)) +
    geom_hline(yintercept = c(11)) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.position = "none",
          plot.title = element_text(color = "Black", size = 14,hjust = 0.75))

  Figure_Mg_free

  output.df = df %>% filter(Sample == "Chelator") %>% select(Metabolites, Conc.Mg, Mg.free)

  output.plot = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, labels = Labels,
                     align = "v")

  df.model$Conc.Mg[which.min(abs(df.model$Mg.free - 2))]

  length(unique(df.model$Conc.Mg))

  output.model = df.model %>% aggregate(Mg.free ~ Conc.Mg, FUN = mean)


  CI95 = c()

  for (i in 1:nrow(output.model)){
    expected.Mgfree = df.model %>%
      filter(Conc.Mg == output.model$Conc.Mg[i])
    CI95[i] = paste(quantile(expected.Mgfree$Mg.free, 0.025), "to", quantile(expected.Mgfree$Mg.free, 0.975))
  }

  output.model$CI95 = CI95

  output.model$Condition = df$Metabolites[1]

  output = list(output.plot,
                output.model,
                output.df,
                output.HQS)

}

df.HQS = read.csv("Figures/Figure_1/HQS_data.csv")

df.AC.model = read.csv("Figures/Figure_1/Modeled_AC_MCM_concentrations.csv")

unique(df.HQS$Metabolites)

df.HQS$Metabolites = factor(df.HQS$Metabolites,
                            levels = c("NTPCM", "WMCM", "Eco80"),
                            labels = c("NTPCM", "WMCM", "Eco80"))

#NTPCM
df.model.NTPCM = df.AC.model %>% select(Mg.T, Mg.free.NTP)
colnames(df.model.NTPCM) = c('Conc.Mg', "Mg.free")

Figure_1BE = analyze.HQS(df.HQS %>% filter(Metabolites == "NTPCM"),
                         df.model.NTPCM,
                         viridis(n =  7)[3],
                         Labels = c("C", "F"),
                         xlimits = c(0, 70))

#WMCM
df.model.WMCM = df.AC.model %>% select(Mg.T, Mg.free.WMCM)
colnames(df.model.WMCM) = c('Conc.Mg', "Mg.free")

Figure_1CF = analyze.HQS(df.HQS %>% filter(Metabolites == "WMCM"),
                         df.model.WMCM,
                         viridis(n =  7)[1],
                         Labels = c("D", "G"),
                         xlimits = c(0, 70))

#Ecoli80
df.model.Ecoli80 = df.AC.model %>% select(Mg.T, Mg.free)
colnames(df.model.Ecoli80) = c('Conc.Mg', "Mg.free")

Figure_1DG = analyze.HQS(df.HQS %>% filter(Metabolites == "Eco80"),
                         df.model.Ecoli80,
                         viridis(n =  7)[6],
                         Labels = c("B", "E"),
                         xlimits = c(1, 80))




Figure_1BCDEFG = plot_grid(Figure_1DG[[1]], Figure_1BE[[1]], Figure_1CF[[1]], nrow = 1)


Figure_1ABCDEFG = plot_grid(Figure_1A, Figure_1BCDEFG, labels = "A", ncol = 1, rel_heights = c(0.5,1.5))

ggsave("Figures/Figure_1/Figure_1.svg", Figure_1ABCDEFG, width = 3.3, height = 3, units = "in", scale = 3)

df.model = bind_rows(Figure_1DG[[2]],Figure_1BE[[2]], Figure_1CF[[2]])

list.files("Figures/Table 2")
write.csv(df.model, "Figures/Table 2/Modeled_results.csv")

df.HQS.fits = bind_rows(Figure_1DG[[4]],Figure_1BE[[4]], Figure_1CF[[4]])
list.files("Figures/SI_Table_3_HQS_fits_in_AC")
write.csv(df.HQS.fits, "Figures/SI_Table_3_HQS_fits_in_AC/HQS_binding_K_fits.csv")


df.final = bind_rows(Figure_1DG[[3]],Figure_1BE[[3]], Figure_1CF[[3]])

write.csv(df.final, "Figures/Table 2/HQS_results.csv")
