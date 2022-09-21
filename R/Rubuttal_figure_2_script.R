
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

df.HQS = read.csv("Figures/Figure_1/HQS_data.csv")

df.AC.model = read.csv("Figures/Figure_1/Modeled_AC_MCM_concentrations.csv")

unique(df.HQS$Metabolites)

#NTPCM
df.model.NTPCM = df.AC.model %>% select(Mg.T, Mg.free.NTP)
colnames(df.model.NTPCM) = c('Conc.Mg', "Mg.free")

#WMCM
df.model.WMCM = df.AC.model %>% select(Mg.T, Mg.free.WMCM)
colnames(df.model.WMCM) = c('Conc.Mg', "Mg.free")

#Ecoli80
df.model.Ecoli80 = df.AC.model %>% select(Mg.T, Mg.free)
colnames(df.model.Ecoli80) = c('Conc.Mg', "Mg.free")


#####Figure 1 B-G####


list.files("Figures/Figure_1")

  df = df.HQS %>% filter(Metabolites == "Eco80")
  df.model = df.model.Ecoli80
  color = viridis(n =  7)[3]
  xlimits = c(0, 75)
  ylimits = c(-1, 20)

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

  fit.form(0.5)

  ####Rebuttal figure 1####

  R_figure_1 = ggplot() +
    geom_polygon(mapping = aes(x = c(0.1, 210, 210, 0.1),
                               y = c(fit.form(c(0.5, 0.5, 3, 3)))),
                 fill = "red",
                 alpha = 0.5) +
    geom_hline(yintercept = c(fit.form(2)), color = "red") +
    geom_point(data = df, mapping = aes(x = Conc.Mg, y = I.norm, shape = Sample, color = Sample)) +
    theme_classic() +
    geom_function(fun = fit.form, color = "dimgrey") +
    ylab("HQS emission") +
    xlab(bquote('Total'~'Mg'^'2+'~'(mM)')) +
    ggtitle(df$Metabolites[1]) +
    scale_color_manual(values = c("dimgrey", "black")) +
    scale_x_continuous(limits = c(0.1, 210), expand = c(0, 0)) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.text = element_text(color = "Black", size = 10),
          legend.title = element_blank(),
          legend.position = c(0.75, 0.76),
          legend.background = element_blank(),
          plot.title = element_text(color = "Black", size = 14,hjust = 0.5))

  R_figure_1

  df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

  errorMg = function(I){
    Imax = coef(fit)[1]
    Imin = coef(fit)[2]
    K = coef(fit)[3]
    SImax = coef(summary(fit))[,2][1]
    SImin = coef(summary(fit))[,2][2]
    SK = coef(summary(fit))[,2][3]
    beta =  (I - Imin)/(Imax - Imin)
    I.minusImin = I - Imin
    Imax.minus.Imin = Imax - Imin
    SI.minusImin = SImin
    SImax.minus.Imin = sqrt(SImax^2 + SImin^2)
    Sbeta = beta*sqrt((SI.minusImin/I.minusImin)^2 + (SImax.minus.Imin/Imax.minus.Imin)^2)
    betaK = beta*K
    SbetaK = betaK*sqrt((Sbeta/beta)^2 + (SK/K)^2)
    Kminus.betaK = K - betaK
    SK.minus.betaK = sqrt(SK^2 + SbetaK^2)
    Mg = beta/(K - beta*K)
    SMg = Mg*sqrt((Sbeta/beta)^2 + (SK.minus.betaK/Kminus.betaK)^2)
  }

  df$freeMgError = errorMg(df$Emission)


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

  df$PMgerror = 100*df$freeMgError/df$Mg.free

  R_figure_2 = ggplot() +
    geom_polygon(mapping = aes(x = c(0.5, 3, 3, 0.5),
                               y = c(0.8, 0.8, 120, 120),
                               fill = "red",
                               alpha = 0.5)) +
    geom_vline(xintercept = c(2), color = "red") +
    geom_hline(yintercept = 5) +
    geom_point(data = df,
                    mapping = aes(x = Mg.free,
                                  y = PMgerror,
                                  color = Sample,
                                  shape = Sample)) +
    theme_classic() +
    ylab(bquote('%error Free'~'Mg'^'2+')) +
    xlab(bquote('Free'~'Mg'^'2+'~'(mM)')) +
    scale_color_manual(values = c("dimgrey", "black")) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 200)) +
    scale_y_continuous(trans = "log10", limits = c(0.8, 120),
                       breaks = c(1, 5, 10, 50, 100),
                       expand = c(0, 0)) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.position = "none",
          plot.title = element_text(color = "Black", size = 14,hjust = 0.75))

  R_figure_2

  R_figure_3 = ggplot() +
    geom_hex(data = df.model %>% filter(Conc.Mg < 200), mapping = aes(x = Conc.Mg, y = Mg.free), bins = 50) +
    geom_polygon(data = NULL, mapping = aes(x = c(0.1, 200, 200, 0.1), y = c(0.5, 0.5, 3, 3)),
                 fill = "red", alpha = 0.25) +
    annotate("segment", y = 2, yend = 2, x = 0.1, xend = 200,
             color = "red") +
    geom_pointrange(data = df.no_edta %>% filter(Sample == "Chelator"),
               mapping = aes(x = Conc.Mg, y = Mg.free,
                             ymin = Mg.free - freeMgError, ymax = Mg.free + freeMgError),
               shape = "triangle") +
    theme_classic() +
    geom_line(data = df.no_edta, mapping = aes(x = Conc.Mg, y = model),  color = "black") +
    scale_fill_viridis(option = "D") +
    #annotate("text", x = 30, y = 0.01, label = paste(round(Mg.total, digits = 2), " mM total Mg2+"), color = "red") +
    #annotate("text", x = 40, y = 5, label = paste(round(free.Mg, digits = 2), " mM free Mg2+"), color = "red") +
    ylab(bquote('Free'~'Mg'^'2+'~'(mM)')) +
    xlab(bquote('Total'~'Mg'^'2+'~'(mM)')) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 210)) +
    scale_y_continuous(lim = ylimits) +
    theme(axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 16),
          axis.title.y = element_text(color = "Black", size = 16),
          legend.position = "none",
          plot.title = element_text(color = "Black", size = 14,hjust = 0.75))

plot_grid(R_figure_1, R_figure_2, R_figure_3, nrow = 1, align =  "h",
          labels = c("A", "B", "C"), label_size = 20)

list.files()

list.files("Reviewer_response_figures")

ggsave("Reviewer_response_figures/R_figure_2.svg", width = 6.5, height = 3, scale = 2)
