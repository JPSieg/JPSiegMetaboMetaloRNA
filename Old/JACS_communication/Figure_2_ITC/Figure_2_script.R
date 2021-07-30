setwd("~/Jacob/Research/Manuscripts/JACS_communication/Figure_2_ITC")

library(tidyverse)
library(viridis)
library(cowplot)
library(MuMIn)

list.files("ITC_data_files")

####Make a file index####

list.files("ITC_data_files")

df = read.csv("ITC_data_index.csv")

####Analyze the data on a loop####

devtools::load_all("~/Jacob/R_packages/MetaboMgITC")

#?read.itc
#?MetaboMgITC

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

####Make Figure A####

Figure_A = list.fits[[9]]$Plot

ggsave("Figure_A.svg", Figure_A)

####Make Figure B####

df.final$ymin = df.final$K - df.final$Std.error.K
df.final$ymax = df.final$K + df.final$Std.error.K
df.final$xmin = df.final$Temp. - df.final$Std.error.Temp.
df.final$xmax = df.final$Temp. + df.final$Std.error.Temp.


df.rect = data.frame(x = c(10, 55, 55, 10), y = c(10, 10, 1/0.002, 1/0.002))

scaleFUN <- function(x){output = round(x, digits = 0)}

Figure_B = ggplot() +
  geom_polygon(data = df.rect, mapping = aes(x = x, y = y), fill = "grey") +
  geom_hline(yintercept = 1/0.002, size = 1.5) +
  stat_smooth(data = df.final, mapping = aes(x = Temp.,
                                             y = K,
                                             color = Metabolite,), method = "lm", se = FALSE) +
  geom_pointrange(data = df.final, mapping = aes(x = Temp., color = Metabolite, y = K, xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)) +
  annotate("text", y = c(600, 403), x = c(27, 27), label = c("Strong", "Weak")) +
  scale_color_manual(values = viridis(10)) +
  scale_y_continuous(trans = "log10", limits = c(10, 15000), expand = c(0, 0), labels = scaleFUN) +
  scale_x_continuous(limits = c(20, 55), expand = c(0, 0)) +
  ylab("K (1/M)") +
  xlab("Temperature (\u00b0C)") + 
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = c(0.9, 0.9))
  


Figure_B

  
  

Figure_B



ggsave("Test.svg", Figure_B)

####Model HQS data####

list.files()

df = read.csv("HQS_titration_model.csv")

K.HQS = (10^-8.43)*(11300)/((10^-7) + 10^-8.43)

log10(K.HQS)

1000/K.HQS

K.HQS = K.HQS/1000 #Make it mM units

HQS.emission = function(Mg, K.HQS = 0.4047952, Em.max = 1, Em.min = 0){
  output = (Em.max - Em.min)*(K.HQS*Mg/(1 + K.HQS*Mg)) + Em.min
}

#No metabolites

df$Emission = HQS.emission(df$Mg.total) + rnorm(length(df$N.injection), 0.0, 0.005)

df$Condition = "None"

#equation for calculating emission with metabolites
HQS.emission.metabolite = function(Mg, Kd, L.total){
  K = 1/Kd #mM
  Xf = L.total #mM
  Mf = Mg #mM
  a = K
  b = K*Xf - K*Mf + 1
  c = -Mf
  Mg.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
  output = Mg.free*K.HQS/(1 + Mg.free*K.HQS)
}

#With strong chelators

df.strong = df

df.strong$Emission = HQS.emission.metabolite(df.strong$Mg.total, 0.2225, 27.41)  + rnorm(length(df$N.injection), 0.0, 0.005)

df.strong$Condition = "Strong"

#With weak chelators

df.weak = df
df.weak$Emission = HQS.emission.metabolite(df.weak$Mg.total, (164-6)*2/6, 164)  + rnorm(length(df$N.injection), 0.0, 0.005)
df.weak$Condition = "Weak"

#With total chelators

df.total = df
df.total$Emission = HQS.emission.metabolite(df.total$Mg.total, (194-27)*2/24, 164)  + rnorm(length(df$N.injection), 0.0, 0.005)
df.total$Condition = "Total"



#plot(df$Mg.total, df$Emission)
#plot(df.strong$Mg.total, df.strong$Emission)
#plot(df.weak$Mg.total, df.weak$Emission)


df.modeled = bind_rows(df, df.strong, df.weak, df.total)

####Make Figure C####

Figure_C = ggplot(df.modeled, aes(x = Mg.total, y = Emission, color = Condition)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  stat_function(fun = HQS.emission, color = "black") +
  ylab("HQS emmission") +
  xlab("[Mg] total (mM)") +
  scale_color_manual(values = viridis(5)) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = c(0.8, 0.3))


ggsave("test.svg", Figure_C)

####Make Figure D####

#none

df$Mg.free = df$Emission/(K.HQS*(1-df$Emission))

list.fit = {}

for (i in 1:10){
  list.fit[[i]] = lm(Mg.free ~  poly(Mg.total, i, raw=TRUE), df)
}

df.model.sel = data.frame(model.sel(list.fit))

best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
best.polynomial.none = list.fit[[best.polynomial.order]]

df$model = predict(list.fit[[best.polynomial.order]])

#weak

df.weak$Mg.free = df.weak$Emission/(K.HQS*(1-df.weak$Emission))

list.fit.weak = {}

for (i in 1:10){
  list.fit.weak[[i]] = lm(Mg.free ~  poly(Mg.total, i, raw=TRUE), df.weak)
}

df.model.sel = data.frame(model.sel(list.fit.weak))

best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
best.polynomial.weak = list.fit.weak[[best.polynomial.order]]

df.weak$model = predict(list.fit.weak[[best.polynomial.order]])

#strong

df.strong$Mg.free = df.strong$Emission/(K.HQS*(1-df.strong$Emission))

list.fit.strong = {}

for (i in 1:10){
  list.fit.strong[[i]] = lm(Mg.free ~  poly(Mg.total, i, raw=TRUE), df.strong)
}

df.model.sel = data.frame(model.sel(list.fit.strong))

best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
best.polynomial.strong = list.fit.strong[[best.polynomial.order]]

df.strong$model = predict(list.fit.strong[[best.polynomial.order]])

#total

df.total$Mg.free = df.total$Emission/(K.HQS*(1-df.total$Emission))

list.fit.total = {}

for (i in 1:10){
  list.fit.total[[i]] = lm(Mg.free ~  poly(Mg.total, i, raw=TRUE), df.total)
}

df.model.sel = data.frame(model.sel(list.fit.total))

best.polynomial.order = as.integer(rownames(df.model.sel)[which.max(df.model.sel$weight)])
best.polynomial.total = list.fit.total[[best.polynomial.order]]

df.total$model = predict(list.fit.total[[best.polynomial.order]])

df.modeled = bind_rows(df, df.weak, df.strong, df.total)

find.Mg.total = function(fit,
                         Mg.total.start = c(1:100)){
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
      Mg.free.error = abs(Mg.free - 2)
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






Mg.total = c(find.Mg.total(best.polynomial.none),
             find.Mg.total(best.polynomial.weak),
             find.Mg.total(best.polynomial.total),
             find.Mg.total(best.polynomial.strong))

options(scipen = 999)

Figure_D = ggplot(df.modeled) +
  geom_abline(slope = 1, intercept = 0, color = "dimgrey", size = 1.0) +
  geom_hline(yintercept = 2, colour = 'black', size = 0.5) +
  geom_vline(xintercept = Mg.total, colour = 'black', size = 0.5) +
  geom_point(mapping = aes(x = Mg.total, y = Mg.free, color = Condition)) +
  geom_line(mapping = aes(x = Mg.total, y = model, color = Condition)) +
  theme_classic() +
  ylab("[Mg] free (mM)") +
  xlab("[Mg] total (mM)") +
  scale_color_manual(values = viridis(5)) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = "none")



ggsave("test.svg", Figure_D)

####Make Figure_E####


df = data.frame("Predicted" = c(2, 8, 23, 29),
                "Experimental" = Mg.total,
                "Condition" = c("None", "Weak", "Strong", "Total"))

Figure_E = ggplot(df, aes(x = Predicted, y = Experimental)) +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1.5) +
  geom_point(size = 3) +
  theme_classic() +
  ylab("HQS titration") +
  xlab("ITC and Calculated") +
  scale_color_manual(values = viridis(5)) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = "none")


#?plot_grid

Figure_C_D_E = plot_grid(Figure_C, Figure_D, Figure_E, nrow = 1, labels = c("C", "D", "E"))

Figure_B_C_D_E = plot_grid(Figure_B, Figure_C_D_E, ncol = 1, labels = c("B"), rel_heights = c(4,3))

Figure_A_B_C_D_E = plot_grid(Figure_A, Figure_B_C_D_E, nrow = 1, labels = c("A"), rel_widths = c(1, 3))

ggsave("Figure_2.svg", Figure_A_B_C_D_E, scale = 3, width = 6.5, height = 3.3, units = "in")

