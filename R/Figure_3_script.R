library(tidyverse)
library(MeltR)
library(viridis)
library(ggbeeswarm)

####Read in data####

vector.files = paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df = lapply(vector.files, read.csv)

####Fit data####

list.fit = {}
list.df.Tms = {}

for (i in 1:length(vector.files)){
  print(i)
  list.fit[[i]] = meltR.F(list.df[[i]],
                          K_error = c(0.65, 0.65))
  list.df.Tms[[i]] = list.fit[[i]]$Tms
  if (list.df[[i]]$Helix[1] == FALSE){
    list.df.Tms[[i]]$Helix = "F"
  }else{
    list.df.Tms[[i]]$Helix = list.df[[i]]$Helix[1]
  }
  list.df.Tms[[i]]$Condition = list.df[[i]]$Condition[1]
}

df.Tms = bind_rows(list.df.Tms)

#####Make Tms plot####

df.Tms$Condition = factor(df.Tms$Condition,
                          levels = c("Monovalent", "NTPCM", "WMCM", "Ecoli80"))

ggplot() +
  scale_color_viridis(name = "[BHQ1] (nM)") +
  facet_wrap(~Helix, nrow = 1) +
  geom_beeswarm(data = df.Tms,
                mapping = aes(x = Condition, y = Tm, color = B)) +
  theme_classic() +
  ylab("Tm (\u00b0C)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 14))

ggsave("Figures/Figure_3/Figure_3.png", width = 7, units = "in")
