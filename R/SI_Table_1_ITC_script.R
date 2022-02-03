
####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(scales)
library(MetaboMgITC2)
options(scipen=10000)

####Load in and fit ITC data####

list.files("Figures/SI_Table_1_ITC")
df = read.csv("Figures/SI_Table_1_ITC/ITC_data_index.csv")

list.fits = {}
list.df.fit = {}

colnames(df)

for (i in 1:length(df$Metabolite)){
  print(paste(df$Metabolite[i], " at ", df$Temperature[i], "C", sep = ""))
  df.Cell = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Cell[i], sep = "/"),
                       df$Syrringe.X[i], df$Syringe.M[i], df$Syringe.C[i],
                       df$Cell.X[i], df$Cell.M[i], df$Cell.C[i]))
  df.blank = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Blank[i], sep = "/"),
                        df$Syrringe.X[i], df$Syringe.M[i], df$Syringe.C[i],
                        df$Cell.X[i], df$Cell.M[i], df$Cell.C[i]))
  list.fits[[i]] = MetaboMgITC2(df.Cell, df.blank,
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
View(df.final)

####Write SI table####

list.files("Figures/SI_Table_1_ITC")
write.csv(df.final, "Figures/SI_Table_1_ITC/Kd_results.csv", row.names = FALSE)
