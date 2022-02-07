library(tidyverse)
library(MeltR)
library(ggrepel)

####Read in data####

vector.files = paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df = lapply(vector.files, read.csv)

####Make a vector of file prefixes to save the data####

vector.prefixes = c()

for (i in 1:length(vector.files)){
  if (list.df[[i]]$Helix[1] == FALSE){
    vector.prefixes[i] = paste("Helix", "F", "in", list.df[[i]]$Condition[1], sep = "_")
  }else{
    vector.prefixes[i] = paste("Helix", list.df[[i]]$Helix[1], "in", list.df[[i]]$Condition[1], sep = "_")
  }
}

####Fit data####

list.fit = {}

list.fit[[1]] = meltR.F(list.df[[1]],
                        K_error_quantile = 0.5,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[1],
                        Save_results = "all")

list.fit[[2]] = meltR.F(list.df[[2]],
                        K_error_quantile = 0.5,
                        K_range = c(50, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[2],
                        Save_results = "all")

list.fit[[3]] = meltR.F(list.df[[3]],
                        K_error_quantile = 0.5,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[3],
                        Save_results = "all")

list.fit[[4]] = meltR.F(list.df[[4]],
                        K_error_quantile = 0.5,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[4],
                        Save_results = "all")

list.fit[[5]] = meltR.F(list.df[[5]],
                        K_error_quantile = 0.25,
                        K_range = c(75, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[5],
                        Save_results = "all")

list.fit[[6]] = meltR.F(list.df[[6]],
                        K_error_quantile = 0.5,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[6],
                        Save_results = "all")

list.fit[[7]] = meltR.F(list.df[[7]],
                        K_error_quantile = 0.5,
                        K_range = c(75, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[7],
                        Save_results = "all")

list.fit[[8]] = meltR.F(list.df[[8]],
                        K_error_quantile = 0.5,
                        K_range = c(75, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[8],
                        Save_results = "all")

list.fit[[9]] = meltR.F(list.df[[9]],
                        K_error_quantile = 0.5,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[9],
                        Save_results = "all")

list.fit[[10]] = meltR.F(list.df[[10]],
                         K_error_quantile = 0.25,
                         K_range = c(100, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[10],
                         Save_results = "all")

list.fit[[11]] = meltR.F(list.df[[11]] %>% filter(!Well %in% c("A7", "A6")),
                         K_error_quantile = 0.25,
                         K_range = c(100, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[11],
                         Save_results = "all")

list.fit[[12]] = meltR.F(list.df[[12]],
                         K_error_quantile = 0.5,
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[12],
                         Save_results = "all")

list.fit[[13]] = meltR.F(list.df[[13]],
                         K_error_quantile = 0.5,
                         K_range = c(75, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[13],
                         Save_results = "all")


list.fit[[14]] = meltR.F(list.df[[14]],
                         K_error_quantile = 0.5,
                         K_range = c(75, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[14],
                         Save_results = "all")

list.fit[[15]] = meltR.F(list.df[[15]],
                         K_error_quantile = 0.5,
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[15],
                         Save_results = "all")

list.fit[[16]] = meltR.F(list.df[[16]],
                         K_error_quantile = 0.5,
                         K_range = c(75, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[16],
                         Save_results = "all")


list.fit[[17]] = meltR.F(list.df[[17]]  %>%
                           filter(!Well %in% c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8")) %>%
                           filter(!Well %in% c("D5", "D6", "D8", "D9", "D10")),
                         K_error_quantile = 0.5,
                         K_range = c(20, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[17],
                         Save_results = "all")

df = list.df[[17]] %>% filter(Reading == 60) %>%
  filter(!Well %in% c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8")) %>%
  filter(!Well %in% c("D5", "D6", "D8", "D9", "D10"))
ggplot(df, aes(x = B, y = Emission, label = Well))+
  geom_point() +
  geom_text_repel()


list.fit[[18]] = meltR.F(list.df[[18]],
                         K_error_quantile = 0.5,
                         K_range = c(200, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[18],
                         Save_results = "all")

list.fit[[19]] = meltR.F(list.df[[19]],
                         K_error_quantile = 0.5,
                         K_range = c(50, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[19],
                         Save_results = "all")

list.fit[[20]] = meltR.F(list.df[[20]],
                         K_error_quantile = 0.5,
                         K_range = c(100, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[20],
                         Save_results = "all")

####Consolidate data from fits####

list.df.Tms = {}
list.df.vh = {}
list.df.Ks = {}

for (i in 1:length(vector.files)){
  print(i)
  list.df.Tms[[i]] = list.fit[[i]]$Tms
  if (list.df[[i]]$Helix[1] == FALSE){
    list.df.Tms[[i]]$Helix = "F"
  }else{
    list.df.Tms[[i]]$Helix = list.df[[i]]$Helix[1]
  }
  list.df.Tms[[i]]$Condition = list.df[[i]]$Condition[1]

  list.df.vh[[i]] = list.fit[[i]]$VantHoff
  if (list.df[[i]]$Helix[1] == FALSE){
    list.df.vh[[i]]$Helix = "F"
  }else{
    list.df.vh[[i]]$Helix = list.df[[i]]$Helix[1]
  }
  list.df.vh[[i]]$Condition = list.df[[i]]$Condition[1]

  list.df.Ks[[i]] = list.fit[[i]]$K
  if (list.df[[i]]$Helix[1] == FALSE){
    list.df.Ks[[i]]$Helix = "F"
  }else{
    list.df.Ks[[i]]$Helix = list.df[[i]]$Helix[1]
  }
  list.df.Ks[[i]]$Condition = list.df[[i]]$Condition[1]
  list.df.Ks[[i]]$K_error = list.df.vh[[i]]$K_error[1]
}

df.Tms = bind_rows(list.df.Tms)
df.vh = bind_rows(list.df.vh)
df.Ks = bind_rows(list.df.Ks)

####Write files for plotting####

write.csv(df.Tms, "Figures/Figure_3/Fits_summary_Tms.csv", row.names = FALSE)
write.csv(df.vh, "Figures/Figure_3/Fits_summary_vh.csv", row.names = FALSE)
write.csv(df.Ks, "Figures/Figure_3/Fits_summary_Ks.csv", row.names = FALSE)
