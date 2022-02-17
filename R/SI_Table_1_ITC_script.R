
####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(scales)
library(MetaboMgITC2)
options(scipen=10000)

####Load in and fit ITC data####

list.files("Figures/SI_Table_1_ITC/ITC_data_files")
df = read.csv("Figures/SI_Table_1_ITC/ITC_data_index.csv")

Metabolite = c()
Temperature = c()
K = c()
Kd = c()
SD.Kd = c()
CI95.Kd = c()
CI95.K = c()
H = c()
CI95.H = c()

list.remove = list(1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   c(1, 26, 29),
                   c(1, 2, 22),
                   c(1, 22, 26, 29))

n.simulations = 1000

list.df.metabolites = {}

for (i in 1:length(df$Metabolite)){
  print(paste(df$Metabolite[i], " at ", df$Temperature[i], "C", sep = ""))

  pb = txtProgressBar(min = 1, max = n.simulations, initial = 1)

  list.df.fit = {}
  for (j in 1:n.simulations) {
    setTxtProgressBar(pb, j)
    df.Cell = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Cell[i], sep = "/"),
                         df$Syrringe.X[i] + rnorm(1, 0, df$SE.Syrringe.X[i]),
                         df$Syringe.M[i] + rnorm(1, 0, df$SE.Syrringe.M[i]),
                         df$Syringe.C[i] + rnorm(1, 0, df$Syringe.C[i]),
                         df$Cell.X[i] + rnorm(1, 0, df$SE.Cell.X[i]),
                         df$Cell.M[i] + rnorm(1, 0, df$SE.Cell.M[i]),
                         df$Cell.C[i]))
    df.blank = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Blank[i], sep = "/"),
                          df$Syrringe.X[i] + rnorm(1, 0, df$SE.Syrringe.X[i]),
                          df$Syringe.M[i] + rnorm(1, 0, df$SE.Syrringe.M[i]),
                          df$Syringe.C[i] + rnorm(1, 0, df$Syringe.C[i]),
                          df$Cell.X[i] + rnorm(1, 0, df$SE.Cell.X[i]),
                          df$Cell.M[i] + rnorm(1, 0, df$SE.Cell.M[i]),
                          df$Cell.C[i]))
    fit = MetaboMgITC2(cell = df.Cell, blank =  df.blank,
                       Fit.start = list(H = df$H.start[i], K = df$K.start[i]),
                      Saturation.threshold = FALSE,
                      Remove.injection = list.remove[[i]])
    list.df.fit[[j]] = fit$Table
  }

  params = bind_rows(list.df.fit)

  Metabolite[i] = df$Metabolite[i]

  Temp = params %>%
    filter(Parameter == "Temp.") %>%
    select(Number)
  Temperature[i] = mean(Temp$Number)

  Ks = params %>%
    filter(Parameter == "K") %>%
    select(Number)
  K[i] = mean(Ks$Number)
  CI95.K[i] = paste(quantile(Ks$Number, 0.025), quantile(Ks$Number, 0.975), sep = " to ")

  Kds = 1000/Ks
  Kd[i] = mean(Kds$Number)
  SD.Kd[i] = sd(Kds$Number)
  CI95.Kd[i] = paste(quantile(Kds$Number, 0.025), quantile(Kds$Number, 0.975), sep = " to ")

  Hs = params %>%
    filter(Parameter == "dH") %>%
    select(Number)
  H[i] = mean(Hs$Number)
  CI95.H[i] = paste(quantile(Hs$Number, 0.025), quantile(Hs$Number, 0.975), sep = " to ")

}

df.final = data.frame(Metabolite,
                      Temperature,
                      K,
                      CI95.K,
                      Kd,
                      SD.Kd,
                      CI95.Kd,
                      H,
                      CI95.H)

####Write SI table####

list.files("Figures/SI_Table_1_ITC")
write.csv(df.final, "Figures/SI_Table_1_ITC/Kd_results.csv", row.names = FALSE)
