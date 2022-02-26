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
list.K.range = {}

list.fit[[1]] = meltR.F(list.df[[1]],
                        K_error_quantile = 0.25,
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[1],
                        Save_results = "all")
list.K.range[[1]] = c(1, 150)

list.fit[[2]] = meltR.F(list.df[[2]],
                        K_error_quantile = 0.75,
                        K_range = c(50, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[2],
                        Save_results = "all")
list.K.range[[2]] = c(50, 750)

list.fit[[3]] = meltR.F(list.df[[3]],
                        K_error_quantile = 1,
                        K_range = c(25, 200),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[3],
                        Save_results = "all")

list.K.range[[3]] = c(25, 200)

list.fit[[4]] = meltR.F(list.df[[4]],
                        K_error_quantile = 0.75,
                        K_range = c(25, 200),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[4],
                        Save_results = "all")

list.K.range[[4]] = c(25, 200)

list.fit[[5]] = meltR.F(list.df[[5]] %>% filter(Temperature < 70),
                        file_path = "Figures/Figure_3/Fit_data/",
                        K_error_quantile = 1,
                        K_range = c(50, 750),
                        file_prefix = vector.prefixes[5],
                        Save_results = "all")

list.K.range[[5]] = c(50, 750)

list.fit[[6]] = meltR.F(list.df[[6]],
                        K_error_quantile = 0.5,
                        K_range = c(10, 150),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[6],
                        Save_results = "all")

list.K.range[[6]] = c(10, 150)

list.fit[[7]] = meltR.F(list.df[[7]],
                        K_error_quantile = 0.75,
                        K_range = c(75, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[7],
                        Save_results = "all")

list.K.range[[7]] = c(75, 750)

list.fit[[8]] = meltR.F(list.df[[8]],
                        K_error_quantile = 0.5,
                        K_range = c(50, 750),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[8],
                        Save_results = "all")

list.K.range[[8]] = c(50, 750)

list.fit[[9]] = meltR.F(list.df[[9]],
                        K_error_quantile = 0.5,
                        K_range = c(15, 150),
                        file_path = "Figures/Figure_3/Fit_data/",
                        file_prefix = vector.prefixes[9],
                        Save_results = "all")

list.K.range[[9]] = c(15, 150)

list.fit[[10]] = meltR.F(list.df[[10]],
                         K_error_quantile = 1,
                         K_range = c(50, 500),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[10],
                         Save_results = "all")

list.K.range[[10]] = c(50, 500)

list.fit[[11]] = meltR.F(list.df[[11]] %>% filter(!Well %in% c("A7", "A6")),
                         K_error_quantile = 0.25,
                         K_range = c(100, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[11],
                         Save_results = "all")

list.K.range[[11]] = c(100, 750)

list.fit[[12]] = meltR.F(list.df[[12]],
                         K_error_quantile = 0.75,
                         K_range = c(20, 200),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[12],
                         Save_results = "all")

list.K.range[[12]] = c(20, 200)

list.fit[[13]] = meltR.F(list.df[[13]],
                         K_error_quantile = 1,
                         K_range = c(75, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[13],
                         Save_results = "all")

list.K.range[[13]] = c(75, 750)

list.fit[[14]] = meltR.F(list.df[[14]],
                         K_error_quantile = 1,
                         K_range = c(75, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[14],
                         Save_results = "all")

list.K.range[[14]] = c(75, 750)

list.fit[[15]] = meltR.F(list.df[[15]],
                         K_error_quantile = 1,
                         K_range = c(25, 300),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[15],
                         Save_results = "all")

list.K.range[[15]] = c(25, 300)

list.fit[[16]] = meltR.F(list.df[[16]],
                         K_error_quantile = 1,
                         K_range = c(25, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[16],
                         Save_results = "all")

list.K.range[[16]] = c(25, 750)

list.fit[[17]] = meltR.F(list.df[[17]]  %>%
                           filter(!Well %in% c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8")) %>%
                           filter(!Well %in% c("D5", "D6", "D8", "D9", "D10")),
                         K_error_quantile = 0.35,
                         K_range = c(20, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[17],
                         Save_results = "all")

list.K.range[[17]] = c(20, 750)

df = list.df[[17]] %>% filter(Reading == 60) %>%
  filter(!Well %in% c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8")) %>%
  filter(!Well %in% c("D5", "D6", "D8", "D9", "D10"))
ggplot(df, aes(x = B, y = Emission, label = Well))+
  geom_point() +
  geom_text_repel()


list.fit[[18]] = meltR.F(list.df[[18]] %>%
                           filter(!Well %in% c("A10", "C10", "C9", "B9", "C8", "A6")),
                         K_error_quantile = 0.3,
                         K_range = c(100, 600),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[18],
                         Save_results = "all")

list.K.range[[18]] = c(100, 600)

df = list.df[[18]] %>% filter(Reading == 1) %>%
  filter(!Well %in% c("A10", "C10", "C9", "B9", "C8", "A6"))
ggplot(df, aes(x = B, y = Emission, label = Well))+
  geom_point() +
  geom_text_repel()

list.fit[[19]] = meltR.F(list.df[[19]],
                         K_error_quantile = 0.75,
                         K_range = c(20, 150),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[19],
                         Save_results = "all")

list.K.range[[19]] = c(20, 150)

list.fit[[20]] = meltR.F(list.df[[20]],
                         K_error_quantile = 0.5,
                         K_range = c(50, 750),
                         file_path = "Figures/Figure_3/Fit_data/",
                         file_prefix = vector.prefixes[20],
                         Save_results = "all")
list.K.range[[20]] = c(50, 750)

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
  list.df.Ks[[i]]$Kd =  (10^9)/list.df.Ks[[i]]$K
  list.df.Ks[[i]]$In_K_range = TRUE
  list.df.Ks[[i]]$In_K_range[which(list.df.Ks[[i]]$Kd < list.K.range[[i]][1])] = FALSE
  list.df.Ks[[i]]$In_K_range[which(list.df.Ks[[i]]$Kd > list.K.range[[i]][2])] = FALSE
}

df.Tms = bind_rows(list.df.Tms)
df.vh = bind_rows(list.df.vh)
df.Ks = bind_rows(list.df.Ks)
df.Ks$SE.lnK = df.Ks$SE.K/df.Ks$K
df.Ks$lnK = log(df.Ks$K)
df.Ks$invT = 1/(df.Ks$Temperature + 273.15)

####Write files for plotting####

write.csv(df.Tms, "Figures/Figure_3/Fits_summary_Tms.csv", row.names = FALSE)
write.csv(df.vh, "Figures/Figure_3/Fits_summary_vh.csv", row.names = FALSE)
write.csv(df.Ks, "Figures/Figure_3/Fits_summary_Ks.csv", row.names = FALSE)

####Plot van't hoff plots####

df.Ks$Condition = factor(df.Ks$Condition,
                         levels = c("Monovalent", "NTPCM", "WMCM", "Ecoli80"))


ggplot(df.Ks %>%
         filter(In_K_range == TRUE) %>%
         filter(SE.lnK <= K_error),
       aes(x = invT, y = lnK,
           ymin = lnK - SE.lnK,
           ymax = lnK +SE.lnK,
           color = Condition,
           shape = Condition)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Helix, scales = "free", nrow = 1) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)])) +
  geom_pointrange() +
  ylab("ln[K (1/M)]") +
  xlab("1/Temperature (K)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        legend.position = c(0.95, 0.3),
        legend.background = element_blank())

list.files("Figures/SI_figure_x_Vant_hoff_plots")

ggsave("Figures/SI_figure_x_Vant_hoff_plots/SI_figure_x_Vant_hoff_plots.png",
       height = 2, scale = 2)
