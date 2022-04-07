library(tidyverse)
library(viridis)
library(cowplot)
library(scales)
library(MuMIn)
library(ggbeeswarm)
library(ggrepel)

####Read in data####

list.files("Figures/SI_Figure_X_HQS_to_measure_weak_binders")

df.HQS = read.csv("Figures/SI_Figure_X_HQS_to_measure_weak_binders/HQS_data.csv")

head(df.HQS)

####Function that analyzes the data for each ligand####

unique(df.HQS$Chelator)

analyze.HQS = function(df = df.HQS,
                       Ligand = "Aspartate"){
  df = df %>% filter(Chelator == Ligand)

  fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
            df %>% filter(Sample == "No chelator"),
            start = list(I.max = 150000, I.min = 0, K = 10))

  fit.form = function(Conc.Mg){
    Emission = (coef(fit)[1] - coef(fit)[2])*(coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg)) + coef(fit)[2]
  }

  Figure_Raw_Em = ggplot(df, aes(x = Conc.Mg, y = Emission, color = Sample)) +
    geom_point() +
    theme_classic() +
    #geom_function(fun = fit.form, color = "black") +
    ylab("HQS emmission") +
    xlab("[Mg] total (mM)") +
    scale_color_manual(values = viridis(5), name = "Chelator") +
    theme(axis.line = element_line(colour = 'black', size = 1.5),
          axis.ticks = element_line(colour = "black", size = 1.5),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 18),
          axis.title.y = element_text(color = "Black", size = 18),
          legend.text = element_text(color = "Black", size = 16),
          legend.title = element_text(color = "Black", size = 16),
          legend.position = c(0.8, 0.3))

  df$I.norm = (df$Emission - coef(fit)[2])/(coef(fit)[1]- coef(fit)[2])

  fit1 = fit
  fit.form = function(Conc.Mg){
    Emission = (coef(fit1)[3]*Conc.Mg/(1 + coef(fit1)[3]*Conc.Mg))
  }

  Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
    geom_point() +
    theme_classic() +
    geom_function(fun = fit.form, color = "black") +
    ylab("HQS emmission") +
    xlab("[Mg] total (mM)") +
    scale_color_manual(values = viridis(5), name = "Ligand") +
    theme(axis.line = element_line(colour = 'black', size = 1.5),
          axis.ticks = element_line(colour = "black", size = 1.5),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 18),
          axis.title.y = element_text(color = "Black", size = 18),
          legend.text = element_text(color = "Black", size = 16),
          legend.title = element_text(color = "Black", size = 16),
          legend.position = c(0.8, 0.3))

  df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

  Mg.free.calc = function(K, Conc.Mg){
    a = K
    b = K*240 - K*Conc.Mg + 1
    c = -Conc.Mg
    x = (-b + sqrt((b^2) - (4*a*c)))/(2*a)
  }

  df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

  Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                     data = df.chelator %>% filter(Conc.Mg >= 0.1) %>% filter(Conc.Mg <= 20),
                     trace = TRUE,
                     algorithm = "port",
                     lower = 0,
                     control = nls.control(warnOnly = TRUE),
                     start = list(K = 0.0000005))

  Mg.free.chelator= function(Conc.Mg){
    a = coef(Fit.chelator)[1]
    b = coef(Fit.chelator)[1]*240 - coef(Fit.chelator)[1]*Conc.Mg + 1
    c = -Conc.Mg
    x = (-b + sqrt((b^2) - (4*a*c)))/(2*a)
  }

  Figure_Mg_free = ggplot(df, mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
    geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
    geom_point() +
    theme_classic() +
    geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
    #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
    scale_color_manual(values = viridis(5)) +
    scale_y_continuous(trans = "log10", labels = comma, limits = c(0.01, 100)) +
    scale_x_continuous(trans = "log10", labels = comma, limits = c(0.01, 100)) +
    ylab("[Mg] free (mM)") +
    xlab("[Mg] total (mM)") +
    theme(axis.line = element_line(colour = 'black', size = 1.5),
          axis.ticks = element_line(colour = "black", size = 1.5),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 18),
          axis.title.y = element_text(color = "Black", size = 18),
          legend.text = element_text(color = "Black", size = 16),
          legend.title = element_text(color = "Black", size = 16),
          legend.position = "none")

  Figure_I.norm_EDTA_check = ggplot(df %>% filter(Conc.Mg < 0.09),
                                    aes(x = Sample, y = I.norm, color = EDTA)) +
    geom_beeswarm() +
    scale_color_manual(values = viridis(5)) +
    ylab("Emission") +
    xlab("") +
    theme_classic() +
    ggtitle("No Mg2+ control")+
    theme(axis.line = element_line(colour = 'black', size = 1.5),
          axis.ticks = element_line(colour = "black", size = 1.5),
          axis.text.x = element_text(color = "Black", size = 16),
          axis.text.y = element_text(color = "Black", size = 16),
          axis.title.x = element_text(color = "Black", size = 18),
          axis.title.y = element_text(color = "Black", size = 18),
          legend.text = element_text(color = "Black", size = 16),
          legend.title = element_text(color = "Black", size = 16),
          legend.position = c(0.9, 0.9))

  Figure_free.Mg_EDTA_check = ggplot(df %>% filter(Conc.Mg < 0.09),
                                     aes(x = Sample, y = Mg.free, color = EDTA)) +
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
          legend.position = c(0.9, 0.9))

  Figure_CD = plot_grid(Figure_I.norm_EDTA_check, Figure_free.Mg_EDTA_check, nrow = 2)

  Figure_ABCD = plot_grid(Figure_Norm_Em, Figure_Mg_free, Figure_CD, nrow = 1)

  ggsave(paste("Glutamine_results_", Temp, "degC.png", sep = ""), Figure_ABCD, width = 12, height = 3, scale = 2)

  output = data.frame(coef(summary(Fit.chelator)))

  output$Temperature = Temp

  output = output
}
