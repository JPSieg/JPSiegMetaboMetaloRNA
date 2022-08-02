
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

list.df.fit = {}
list.plot = {}

for (i in 1:length(df$Metabolite)){
  print(paste(df$Metabolite[i], " at ", df$Temperature[i], "C", sep = ""))

  pb = txtProgressBar(min = 1, max = n.simulations, initial = 1)

  df.Cell = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Cell[i], sep = "/"),
                         df$Syrringe.X[i],
                         df$Syringe.M[i],
                         df$Syringe.C[i],
                         df$Cell.X[i],
                         df$Cell.M[i],
                         df$Cell.C[i]))
    df.blank = read.itc(c(paste("Figures/SI_Table_1_ITC/ITC_data_files", df$Blank[i], sep = "/"),
                          df$Syrringe.X[i],
                          df$Syringe.M[i],
                          df$Syringe.C[i],
                          df$Cell.X[i],
                          df$Cell.M[i],
                          df$Cell.C[i]))
    fit = MetaboMgITC2(cell = df.Cell, blank =  df.blank,
                       Fit.start = list(H = df$H.start[i], K = df$K.start[i]),
                      Saturation.threshold = FALSE,
                      Remove.injection = list.remove[[i]])
    list.df.fit[[i]] = fit$Table
    list.df.fit[[i]]$Metabolite = df$Metabolite[i]

    list.plot[[i]] = fit$Plot

}

df.final = bind_rows(list.df.fit)

####Write SI table####

list.files("Figures/SI_Table_1_ITC")
write.csv(df.final, "Figures/SI_Table_1_ITC/Kd_results.csv", row.names = FALSE)

####Write SI figure####

list.files("Figures/SI_Figure_X_ITC")

length(list.plot)

vector.labels = paste(c("(A)", "(B)", "(C)", "(D)", "(E)",
                        "(F)", "(G)", "(H)", "(I)", "(J)"),
                      df$Metabolite[c(4,5,6,7,8,2,9,10,3,1)])

Figure_x = plot_grid(list.plot[[4]],
          list.plot[[5]],
          list.plot[[6]],
          list.plot[[7]],
          list.plot[[8]],
          list.plot[[2]],
          list.plot[[9]],
          list.plot[[10]],
          list.plot[[3]],
          list.plot[[1]], ncol = 5, hjust = 0.075, label_size = 20,
          labels = c("(A)", "(B)", "(C)", "(D)", "(E)",
                     "(F)", "(G)", "(H)", "(I)", "(J)"))

ggsave("Figures/SI_Figure_X_ITC/SI_Figure_X_ITC.svg",
       Figure_x, width = 6, height = 4, scale = 4)
ggsave("Figures/SI_Figure_X_ITC/SI_Figure_X_ITC.png",
       Figure_x, width = 6, height = 4, scale = 4)

