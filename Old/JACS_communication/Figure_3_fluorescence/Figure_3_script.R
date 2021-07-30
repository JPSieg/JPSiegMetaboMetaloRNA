setwd("~/Jacob/Research/Manuscripts/JACS_communication/Figure_3_fluorescence")

library(tidyverse)
library(viridis)
library(cowplot)
library(MeltR)
library(ggbeeswarm)

####Load data####

list.files()

df.index = read.csv("Data_index.csv")

file.path.vector = paste("Fluorescence_data", df.index$Experiment, sep = "/")

list.data = lapply(file.path.vector, read.csv)

####Fit data####

list.fits <- {}

for (i in 1:length(list.data)){
  tryCatch({
    print(i)
    list.fits[[i]] <- meltR.F(list.data[[i]],
                         K_error = c(0.6, 0.6),
                         low_K = 0.01)
  }, error = function(e){
    print(paste("fit", i, "failed"))
  })
}

list.fits[[11]] <- meltR.F(list.data[[11]],
                      K_error = c(0.75, 0.75),
                      low_K = 0.01)

list.fits[[12]] <- meltR.F(list.data[[12]],
                           K_error = c(1, 1),
                           low_K = 0.01)

####Compile Tms####

list.Tms <- {}

for (i in 1:length(list.fits)){
  list.Tms[[i]] <- list.fits[[i]]$Tms
  list.Tms[[i]]$Helix <- df.index$Helix[i]
  list.Tms[[i]]$Condition <- df.index$Condition[i]
}

df.Tms <- bind_rows(list.Tms)

####Make Figure 3A####

df.Tms$Condition = factor(df.Tms$Condition,
                          levels = c("1 M NaCL", "Prokaryotic salts", "Strong MCM", "Weak MCM", "Total MCM", "EDTA","free Mg"))

Figure_A = ggplot(df.Tms, aes(x = Condition, y = Tm, color = Helix)) +
  geom_beeswarm() +
  scale_colour_manual(values = viridis(10)) +
  theme_classic() +
  scale_x_discrete(labels=c("1 M NaCL", "Prokaryotic salts", "Strong MCM", "Weak MCM", "Total MCM", "EDTA","free Mg")) +
  ylab("Tm (\u00b0C)") +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 16, angle = 15, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 18),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16))

ggsave("Figure_3.svg", Figure_A, width = 3.3, height = 3.3, units = "in", scale = 2)
