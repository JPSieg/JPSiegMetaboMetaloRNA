library(tidyverse)
library(MeltR)
library(ggrepel)

####Read in data####

list.files("Figures/SI_figure_x_1MNaCl_energies")

df = read.csv("Figures/SI_figure_x_1MNaCl_energies/10_js5003_Helix_J.csv")


####Fit data####

meltR.F(df,
      Kd_error_quantile = 0.35,
      file_path = "Figures/SI_figure_x_1MNaCl_energies/",
      file_prefix = "Helix_J_1MNaCl",
      Save_results = "all")

