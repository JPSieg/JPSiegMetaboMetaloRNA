setwd("~/Jacob/Research/Manuscripts/JPSiegMetaboMetaloRNA/Figures/Figure_2")

library(tidyverse)
library(viridis)
library(cowplot)
library(MuMIn)

####Load in data and split into a list of data.frames####

df = read.csv("Top_15_E.coli_Metabolites_conc_Kds.csv")

head(df)

list.df.T = {}

for (i in 1:length(unique(df$Temperature))){
  list.df.T[[i]] = df %>% filter(Temperature == unique(df$Temperature)[i])
}

names(list.df.T) = unique(df$Temperature)

####Calculate 2 mM free Mg at 37C####

list.df.T$`37`$MCM = list.df.T$`37`$Concentration*2/(list.df.T$`37`$Kd + 2)

list.df.T$`37`$MCM[length(list.df.T$`37`$MCM)] = 2

Mg.Total = sum(list.df.T$`37`$MCM, na.rm = TRUE)

####Numerically determine how free Mg changes with temperature####

Numeric.free.Mg.solver = function(df){
  go = TRUE
  First.recursion = TRUE

  while(go){
    if (First.recursion){
      Mg.free = c(1:Mg.Total)
      First.recursion = FALSE
    }

    error.Mg.T = c()

    for (i in 1:length(Mg.free)){
      MCM = df$Concentration*Mg.free[i]/(df$Kd + Mg.free[i])
      error.Mg.T[i] = abs(Mg.Total - Mg.free[i] - sum(MCM, na.rm = TRUE))
    }

    center = Mg.free[which.min(error.Mg.T)]

    error = error.Mg.T[which.min(error.Mg.T)]

    Mg.free = seq(center - error, center + error, length.out = 50)

    if (error <= 0.00001){go = FALSE}
  }
  print(center)
  output = center
}

Free.Mg = lapply(list.df.T, Numeric.free.Mg.solver)

####Consolidate data####

for (i in 1:length(list.df.T)){
  list.df.T[[i]]$MCM = list.df.T[[i]]$Concentration*Free.Mg[[i]]/(list.df.T[[i]]$Kd + Free.Mg[[i]])
  list.df.T[[i]]$MCM[length(list.df.T[[i]]$MCM)] = Free.Mg[[i]]
}

df = bind_rows(list.df.T)

Strength = c()

for(i in 1:length(df$Metabolites)){
  if (is.na(df$Kd[i]) == FALSE){
    if(df$Kd[i] >= 2){Strength[i] = "Weak"}
    if(df$Kd[i] < 2){Strength[i] = "Strong"}
  }else{
    if(df$Metabolites[i] == "Free Mg"){Strength[i] = "Free Mg"}
  }
}

df$Strength = Strength

####Make Figure 2A####

Figure_2A = ggplot(df, aes(x = factor(Temperature), y = MCM, fill = Strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  theme_minimal()+
  ylab("[Magnesium] (M)") +
  xlab("Temperature (\u00b0C)") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  theme(axis.line.y = element_line(colour = 'black', size = 1.5),
        axis.line.x = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.ticks.x =element_line(colour = "black", size = 1.5),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        legend.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16))

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

####Make Figure D####

Figure_D = ggplot(df.modeled, aes(x = Mg.total, y = Emission, color = Condition)) +
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

####Make Figure E####

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

Figure_E = ggplot(df.modeled) +
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


####Make Figure_F####


df = data.frame("Predicted" = c(2, 8, 23, 29),
                "Experimental" = Mg.total,
                "Condition" = c("None", "Weak", "Strong", "Total"),
                "Temperature" = c(37))

Figure_F = ggplot(df, aes(x = Predicted, y = Experimental, shape = Condition, color = Temperature)) +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1.5) +
  geom_point(size = 4) +
  theme_classic() +
  ylab("Mg total HQS titration (mM)") +
  xlab("Mg total ITC and Calculated (mM)") +
  scale_color_viridis() +
  coord_fixed() +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 18),
        axis.title.y = element_text(color = "Black", size = 18),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16),
        legend.position = c(0.15, 0.71))

Figure_DEF = plot_grid(Figure_D, Figure_E, Figure_F, labels = c("D", "E", "F"), label_size = 20, nrow = 1)

Figure_ADEF =plot_grid(Figure_2A, Figure_DEF, labels = c("A"), label_size = 20, nrow = 2)

ggsave("Figure_2.svg", Figure_ADEF, scale = 5, width = 3.3, height = 2, units = "in", bg = "white")
