library(R2easyR)
library(tidyverse)

####Gua aptamer####

vector.files = paste("Figures/SI_figure_x_R2easyR/Gua_data_files", list.files("Figures/SI_figure_x_R2easyR/Gua_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)
df = bind_rows(list.df)
hist(df$Reactivity, breaks =20)
abline(v = c(100, 1000))

setwd("Figures/SI_figure_x_R2easyR")

palettes = r2easyR.palettes()

vector.files

Condition = c("Published", "25mMFree", "2mMFree", "Eco80", "NTPCM", "WMCM")

table(list.df[[1]]$Dotbracket)

for (i in 2:length(list.df)){
  list.df[[i]]$Dotbracket[c(13, 28, 29, 37, 43, 44, 51, 52)] = "."

  list.df[[i]]$Reactivity[list.df[[i]]$Reactivity <= 0] = NA
  list.df[[i]] = r2easyR.color(list.df[[i]],
                               palettes$Reds.c,
                               manual.scale = c(100, 1000))

  prefix = gsub("Figures/SI_figure_x_R2easyR/Gua_data_files", "R2R", strsplit(vector.files[i], split = "[.]")[[1]][1])

  r2easyR.write(prefix, list.df[[i]], RNA_name = Condition[i], colors = "circles")
}

r2easyR.write("Figures/Figure_2/R2R/No_reactivity", list.df[[1]], RNA_name = "No_reactivity", colors = "NA")



i = 1

list.df[[i]]$Reactivity[list.df[[i]]$Reactivity <= 0] = NA
list.df[[i]] = r2easyR.color(list.df[[i]] %>% filter(N <= 73),
                             palettes$Reds.c,
                             manual.scale = c(4000, 40000))

prefix = gsub("Figures/SI_figure_x_R2easyR/Gua_data_files", "R2R", strsplit(vector.files[i], split = "[.]")[[1]][1])

r2easyR.write(prefix, list.df[[i]], RNA_name = Condition[i], colors = "circles")
