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

####Analyze Glutamate data ####

unique(df.HQS$Chelator)

Ligand = "Glutamate"
df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

Figure_Raw_Em = ggplot(df, aes(x = Conc.Mg, y = Emission, color = Sample)) +
    geom_point() +
    theme_classic() +
    #geom_function(fun = fit.form, color = "black") +
  ggtitle("L-Glutamate") +
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

fit.form = function(Conc.Mg){
    Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
    geom_point() +
    theme_classic() +
    geom_function(fun = fit.form, color = "black") +
  ggtitle("L-Glutamate") +
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
          legend.position = c(0.5, 0.3),
          plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
    a = 1
    b = 240 + Conc.Mg + (1/K)
    c = 240*Conc.Mg
    x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                     data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 50),
                     trace = TRUE,
                     algorithm = "port",
                     lower = 0,
                     control = nls.control(warnOnly = TRUE),
                     start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
    a = 1
    b = 240 + Conc.Mg + (1/coef(Fit.chelator)[1])
    c = 240*Conc.Mg
    x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
  }

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
    geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
    geom_point() +
    theme_classic() +
    geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
    #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
    scale_color_manual(values = viridis(5)) +
    scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
    scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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


Figure_Glutamate = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_glutamate.csv",
          row.names = FALSE)

df.result = data.frame(coef(summary(Fit.chelator)))
df.result = bind_rows(df.result,
                      data.frame("Estimate" = (1/df.result$Estimate),
                                 "Std..Error" = 1/df.result$Estimate*(df.result$Std..Error/df.result$Estimate),
                                 "t.value" = NA,
                                 "Pr...t.." = NA))

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_glutamate.csv",
          row.names = FALSE)

####Analyze Glutathione data ####

unique(df.HQS$Chelator)

Ligand = "Glutathione"
df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

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
        legend.position = c(0.8, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$I.norm = (df$Emission - coef(fit)[2])/(coef(fit)[1]- coef(fit)[2])

fit.form = function(Conc.Mg){
  Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
  geom_point() +
  theme_classic() +
  geom_function(fun = fit.form, color = "black") +
  ggtitle("Glutathione") +
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
        legend.position = c(0.5, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
  a = 1
  b = 193.8 + Conc.Mg + (1/K)
  c = 193.8*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                   data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 60),
                   trace = TRUE,
                   algorithm = "port",
                   control = nls.control(warnOnly = TRUE),
                   start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
  a = 1
  b = 193.8 + Conc.Mg + (1/coef(Fit.chelator))
  c = 193.8*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_point() +
  theme_classic() +
  geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
  #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
  scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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


Figure_Glutathione = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_glutathione.csv",
          row.names = FALSE)

df.result = data.frame(coef(summary(Fit.chelator)))
df.result = bind_rows(df.result,
                      data.frame("Estimate" = (1/df.result$Estimate),
                                 "Std..Error" = 1/df.result$Estimate*(df.result$Std..Error/df.result$Estimate),
                                 "t.value" = NA,
                                 "Pr...t.." = NA))

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_glutathione.csv",
          row.names = FALSE)

####Analyze Aspartate data ####

unique(df.HQS$Chelator)

Ligand = "Aspartate"

df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

Figure_Raw_Em = ggplot(df, aes(x = Conc.Mg, y = Emission, color = Sample)) +
  geom_point() +
  theme_classic() +
  ggtitle("L-Aspartate") +
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

fit.form = function(Conc.Mg){
  Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
  geom_point() +
  theme_classic() +
  geom_function(fun = fit.form, color = "black") +
  ggtitle("L-Aspartate") +
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
        legend.position = c(0.5, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/K)
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                   data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 20),
                   trace = TRUE,
                   algorithm = "port",
                   lower = 0,
                   control = nls.control(warnOnly = TRUE),
                   start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/coef(Fit.chelator)[1])
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_point() +
  theme_classic() +
  geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
  #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
  scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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


Figure_Aspartate = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_aspartate.csv",
          row.names = FALSE)

df.result = data.frame(coef(summary(Fit.chelator)))
df.result = bind_rows(df.result,
                      data.frame("Estimate" = (1/df.result$Estimate),
                                 "Std..Error" = 1/df.result$Estimate*(df.result$Std..Error/df.result$Estimate),
                                 "t.value" = NA,
                                 "Pr...t.." = NA))

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_aspartate.csv",
          row.names = FALSE)

####Analyze Valine data ####

unique(df.HQS$Chelator)

Ligand = "Valine"

df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 0.4),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

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

fit.form = function(Conc.Mg){
  Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
  geom_point() +
  theme_classic() +
  geom_function(fun = fit.form, color = "black") +
  ylab("HQS emmission") +
  xlab("[Mg] total (mM)") +
  ggtitle("L-Valine") +
  scale_color_manual(values = viridis(5), name = "Ligand") +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = c(0.5, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/K)
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                   data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 20),
                   trace = TRUE,
                   algorithm = "port",
                   control = nls.control(warnOnly = TRUE),
                   start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/coef(Fit.chelator)[1])
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_point() +
  theme_classic() +
  geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
  #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
  scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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

Figure_Valine = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_valine.csv",
          row.names = FALSE)

df.result = data.frame(coef(summary(Fit.chelator)))
df.result = bind_rows(df.result,
                      data.frame("Estimate" = (1/df.result$Estimate),
                                 "Std..Error" = 1/df.result$Estimate*(df.result$Std..Error/df.result$Estimate),
                                 "t.value" = NA,
                                 "Pr...t.." = NA))

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_valine.csv",
          row.names = FALSE)


####Analyze Glutamine data####

unique(df.HQS$Chelator)

Ligand = "Glutamine"

df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

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

fit.form = function(Conc.Mg){
  Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
  geom_point() +
  theme_classic() +
  geom_function(fun = fit.form, color = "black") +
  ggtitle("L-Glutamine") +
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
        legend.position = c(0.5, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/K)
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                   data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 50)%>% filter(Conc.Mg >= 0.1),
                   trace = TRUE,
                   algorithm = "port",
                   control = nls.control(warnOnly = TRUE),
                   start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
  a = 1
  b = 240 + Conc.Mg + (1/coef(Fit.chelator)[1])
  c = 240*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_point() +
  theme_classic() +
  geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
  #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
  scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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


Figure_Glutamine = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_glutamine.csv",
          row.names = FALSE)

df.result = data.frame(NA)

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_glutamine.csv",
          row.names = FALSE)

####Analyze pyruvate data####

unique(df.HQS$Chelator)

Ligand = "Pyruvate"

df = df.HQS %>% filter(Chelator == Ligand)

fit = nls(Emission ~ (I.max - I.min)*(K*Conc.Mg/(1 + K*Conc.Mg)) + I.min,
          df %>% filter(Sample == "No chelator") %>% filter(EDTA == "EDTA = 0 mM"),
          start = list(I.max = 150000, I.min = 0, K = 10),
          algorithm = "port",
          lower = c(0, 0, 0),
          control = nls.control(warnOnly = TRUE))

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

fit.form = function(Conc.Mg){
  Emission = (coef(fit)[3]*Conc.Mg/(1 + coef(fit)[3]*Conc.Mg))
}

Figure_Norm_Em = ggplot(df, aes(x = Conc.Mg, y = I.norm, color =  Sample)) +
  geom_point() +
  theme_classic() +
  geom_function(fun = fit.form, color = "black") +
  ggtitle("Pyruvic acid") +
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
        legend.position = c(0.5, 0.3),
        plot.title = element_text(hjust = 0.5, color = "Black", size = 16))

df$Mg.free = df$I.norm/(coef(fit)[3]*(1 - df$I.norm))

Mg.free.calc = function(K, Conc.Mg){
  a = 1
  b = 5 + Conc.Mg + (1/K)
  c = 5*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

df.chelator =  df %>% filter(Sample == "Chelator") %>% filter(EDTA != "EDTA = 80 mM")

Fit.chelator = nls(Mg.free ~ Mg.free.calc(K, Conc.Mg),
                   data = df.chelator %>% filter(EDTA == "EDTA = 0 mM") %>% filter(Conc.Mg <= 3)%>% filter(Conc.Mg >= 0.1),
                   trace = TRUE,
                   algorithm = "port",
                   control = nls.control(warnOnly = TRUE),
                   start = list(K = 0.000001))

Mg.free.chelator = function(Conc.Mg){
  a = 1
  b = 5 + Conc.Mg + (1/coef(Fit.chelator)[1])
  c = 5*Conc.Mg
  x = Conc.Mg - (-b + sqrt((b^2) + (4*a*c)))/(2*a)
}

Figure_Mg_free = ggplot(df %>% filter(EDTA == "EDTA = 0 mM"), mapping = aes(x = Conc.Mg, y = Mg.free, color = Sample)) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_point() +
  theme_classic() +
  geom_function(fun = Mg.free.chelator, color = viridis(5)[1]) +
  #geom_function(fun = Mg.free.Fru, color = viridis(5)[1]) +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
  scale_x_continuous(trans = "log10", labels = comma, limits = c(0.05, 100)) +
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


Figure_Pyruvate = plot_grid(Figure_Norm_Em, Figure_Mg_free, ncol = 1, align = "v")

write.csv(data.frame(coef(summary(fit))),
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/No_chelator_fit_pyruvate.csv",
          row.names = FALSE)

df.result = data.frame(coef(summary(Fit.chelator)))
df.result = bind_rows(df.result,
                      data.frame("Estimate" = 1/df.result$Estimate,
                                 "Std..Error" = (1/df.result$Estimate)*(df.result$Std..Error/df.result$Estimate),
                                 "t.value" = NA,
                                 "Pr...t.." = NA))

write.csv(df.result,
          "Figures/SI_Figure_X_HQS_to_measure_weak_binders/Kapp_fit_pyruvate.csv",
          row.names = FALSE)

####Make plot####

plot_grid(Figure_Glutamate, Figure_Glutathione,
          Figure_Aspartate, Figure_Valine,
          Figure_Glutamine, Figure_Pyruvate,
          labels = c("A", "B", "C", "D", "E", "F"),
          label_size = 20)

ggsave("Figures/SI_Figure_X_HQS_to_measure_weak_binders/SI_Figure_X_HQS_binding.svg",
       width = 6, scale = 3)

ggsave("Figures/SI_Figure_X_HQS_to_measure_weak_binders/SI_Figure_X_HQS_binding.png",
       width = 6, scale = 3)
