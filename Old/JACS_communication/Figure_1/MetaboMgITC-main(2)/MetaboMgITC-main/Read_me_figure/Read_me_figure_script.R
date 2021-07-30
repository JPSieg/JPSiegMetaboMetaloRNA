setwd("~/Jacob/R_packages/MetaboMgITC")

devtools::load_all()

####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)

####Load in data####

#list.files("Read_me_figure")

df <- read.csv("Read_me_figure/Table_1_final.csv")

####Make E.coli figures####

#head(df)

E.coli <- subset(df, Species == "E.coli")

#colnames(df)

#Calculate the Sum of the metabolites

a <- E.coli$Concentration[1]

for (i in 2:length(E.coli$Concentration)){
  a <- c(a, a[i-1] + E.coli$Concentration[i])
}

E.coli$Metabolites.sum <- a

E.coli$Metabolites.sum

#colnames(E.coli)

#plot(E.coli$Concentration, E.coli$MCM)
#abline(a = 0, b = 1)

#Calculate the Sum of the MCM

a <- E.coli$MCM[1]

for (i in 2:length(E.coli$Concentration)){
  if (is.na(E.coli$MCM[i])){
    a <- c(a, a[i-1] + 0)
  }
  if (is.na(E.coli$MCM[i]) == FALSE){
    a <- c(a, a[i-1] + E.coli$MCM[i])
  }
}

E.coli$MCM.sum <- a

####Identify Mg binding strength####

#range(E.coli$Kd.app, na.rm = TRUE)

a <- c()

for (i in 1:length(E.coli$Kd.app)){
  print(i)
  if (is.na(E.coli$Kd.app[i])){a[i] <- "No chelation"}
  if (is.na(E.coli$Kd.app[i]) == FALSE){
    if (E.coli$Kd.app[i] >= 2){a[i] <- "weak"}
    if (E.coli$Kd.app[i] < 2){a[i] <- "strong"}
  }
}

E.coli$Mg.binding.strength <- a

####Figure A####

#colnames(E.coli)

Figure_A <- ggplot(data = E.coli) +
  geom_bar(mapping = aes(x = Metabolites, y = Concentration),
           stat = "identity", width = 1, fill = "#440154FF") +
  geom_bar(mapping = aes(x = Metabolites, y = MCM),
           stat = "identity", width = 1, fill = "#35B779FF") +
  scale_x_discrete(limits = E.coli$Metabolites) +
  geom_point(aes(x= Metabolites, y= Metabolites.sum), stat="identity", color = "#440154FF")+
  geom_line(aes(x= Metabolites, y= Metabolites.sum, group = 1),
            stat="identity", color = "#440154FF", size = 1)+
  geom_point(aes(x= Metabolites, y= MCM.sum), stat="identity", color = "#35B779FF")+
  geom_line(aes(x= Metabolites, y= MCM.sum, group = 1),
            stat="identity", color = "#35B779FF", size = 1)+
  geom_hline(yintercept = sum(E.coli$Concentration, na.rm = TRUE), color =  "#440154FF", size = 1) +
  annotate("text", x = 15, y = 245, label = paste("Estimated total metabolites = ", round(sum(E.coli$Concentration, na.rm = TRUE))," mM", sep =""), size = 4, color = "#440154FF") +
  geom_hline(yintercept = sum(E.coli$MCM, na.rm = TRUE), color =  "#35B779FF", size = 1) +
  annotate("text", x = 25, y = 60, label = paste("Estimated total chelated Mg2+ = ", round(sum(E.coli$MCM, na.rm = TRUE))," mM", sep =""), size = 4, color = "#35B779FF") +
  xlim(E.coli$Metabolites[1:40]) +
  ylim(0,250) +
  ylab("Concentration (mM)")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 12),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16))

#ggsave("Test.png", Figure_C)

####Make Figure B####

#colnames(E.coli)

E.coli.1 = E.coli %>% filter(Kd.app < 2)
E.coli.2 = E.coli %>% filter(Kd.app >= 2)
#View(E.coli.2)

Figure_B <- ggplot(E.coli, aes(x = "", y = Concentration, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  annotate("text", x = 1.75, y = 30, label = E.coli.2$Metabolites[1], size=3) +
  annotate("text", x = 1.63, y = 95, label = E.coli.2$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 125, label = E.coli.2$Metabolites[3], size=3) +
  annotate("text", x = 1.65, y = 205, label = E.coli.1$Metabolites[1], size=3) +
  annotate("text", x = 1.65, y = 215, label = E.coli.1$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 225, label = E.coli.1$Metabolites[3], size=3) +
  annotate("text", x = 1.8, y = 0, label = paste("Total Metabolite =", round(sum(E.coli$Concentration)), "mM"), size=6) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_C)

####Make Figure C####

Figure_C <- ggplot(E.coli, aes(x = "", y = MCM, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(4,6)]) +
  annotate("text", x = 1.67, y = 3, label = E.coli.2$Metabolites[1], size=3) +
  annotate("text", x = 1.67, y = 6, label = E.coli.2$Metabolites[2], size=3) +
  annotate("text", x = 1.72, y = 8, label = E.coli.2$Metabolites[3], size=3) +
  annotate("text", x = 1.72, y = 11, label = E.coli.2$Metabolites[4], size=3) +
  annotate("text", x = 1.75, y = 13, label = E.coli.2$Metabolites[5], size=3) +
  annotate("text", x = 1.67, y = 20, label = E.coli.1$Metabolites[1], size=3) +
  annotate("text", x = 1.67, y = 30, label = E.coli.1$Metabolites[2], size=3) +
  annotate("text", x = 1.67, y = 37, label = E.coli.1$Metabolites[3], size=3) +
  annotate("text", x = 1.67, y = 42, label = E.coli.1$Metabolites[4], size=3) +
  annotate("text", x = 1.8, y = 0, label = paste("Total MCM =", round(sum(E.coli$MCM, na.rm = TRUE)), "mM"), size=6) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_C)

####Make Yeast figures####

#head(df)

Yeast <- subset(df, Species == "Yeast")

#colnames(Yeast)

#head(Yeast)

#Calculate the Sum of the metabolites

a <- Yeast$Concentration[1]

for (i in 2:length(Yeast$Concentration)){
  a <- c(a, a[i-1] + Yeast$Concentration[i])
}

Yeast$Metabolites.sum <- a

Yeast$Metabolites.sum

#colnames(Yeast)

#plot(Yeast$Concentration, Yeast$MCM)
#abline(a = 0, b = 1)

#Calculate the Sum of the MCM

a <- Yeast$MCM[1]

for (i in 2:length(Yeast$Concentration)){
  if (is.na(Yeast$MCM[i])){
    a <- c(a, a[i-1] + 0)
  }
  if (is.na(Yeast$MCM[i]) == FALSE){
    a <- c(a, a[i-1] + Yeast$MCM[i])
  }
}

Yeast$MCM.sum <- a

####Identify Mg binding strength####

#range(Yeast$Kd.app, na.rm = TRUE)

a <- c()

for (i in 1:length(Yeast$Kd.app)){
  print(i)
  if (is.na(Yeast$Kd.app[i])){a[i] <- "No chelation"}
  if (is.na(Yeast$Kd.app[i]) == FALSE){
    if (Yeast$Kd.app[i] >= 2){a[i] <- "weak"}
    if (Yeast$Kd.app[i] < 2){a[i] <- "strong"}
  }
}

Yeast$Mg.binding.strength <- a


####Figure D####

#colnames(Yeast)

Figure_D <- ggplot(data = Yeast) +
  geom_bar(mapping = aes(x = Metabolites, y = Concentration),
           stat = "identity", width = 1, fill = "#440154FF") +
  geom_bar(mapping = aes(x = Metabolites, y = MCM),
           stat = "identity", width = 1, fill = "#35B779FF") +
  scale_x_discrete(limits = Yeast$Metabolites) +
  geom_point(aes(x= Metabolites, y= Metabolites.sum), stat="identity", color = "#440154FF")+
  geom_line(aes(x= Metabolites, y= Metabolites.sum, group = 1),
            stat="identity", color = "#440154FF", size = 1)+
  geom_point(aes(x= Metabolites, y= MCM.sum), stat="identity", color = "#35B779FF")+
  geom_line(aes(x= Metabolites, y= MCM.sum, group = 1),
            stat="identity", color = "#35B779FF", size = 1)+
  geom_hline(yintercept = sum(Yeast$Concentration, na.rm = TRUE), color =  "#440154FF", size = 1) +
  annotate("text", x = 15, y = 250, label = paste("Estimated total metabolites = ", round(sum(Yeast$Concentration, na.rm = TRUE))," mM", sep =""), size = 4, color = "#440154FF") +
  geom_hline(yintercept = sum(Yeast$MCM, na.rm = TRUE), color =  "#35B779FF", size = 1) +
  annotate("text", x = 25, y = 57, label = paste("Estimated total chelated Mg2+ = ", round(sum(Yeast$MCM, na.rm = TRUE))," mM", sep =""), size = 4, color = "#35B779FF") +
  xlim(Yeast$Metabolites[1:40]) +
  ylim(0,250) +
  ylab("Concentration (mM)")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 12),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16))

#ggsave("Test.png", Figure_D)


####Make Figure E####

#colnames(Yeast)

Yeast.1 = Yeast %>% filter(Kd.app < 2)
Yeast.2 = Yeast %>% filter(Kd.app >= 2)
#View(Yeast.2)

Figure_E <- ggplot(Yeast, aes(x = "", y = Concentration, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  annotate("text", x = 1.75, y = 30, label = Yeast.2$Metabolites[1], size=3) +
  annotate("text", x = 1.75, y = 50, label = Yeast.2$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 80, label = Yeast.2$Metabolites[3], size=3) +
  annotate("text", x = 1.65, y = 115, label = Yeast.2$Metabolites[4], size=3) +
  annotate("text", x = 1.65, y = 130, label = Yeast.2$Metabolites[5], size=3) +
  annotate("text", x = 1.8, y = 0, label = paste("Total Metabolite =", round(sum(Yeast$Concentration)), "mM"), size=6) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_E)

####Make Figure F####

Figure_F <- ggplot(Yeast, aes(x = "", y = MCM, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(4,6)]) +
  annotate("text", x = 1.67, y = 0.3, label = Yeast.2$Metabolites[1], size=3) +
  annotate("text", x = 1.7, y = 1, label = Yeast.2$Metabolites[10], size=3) +
  annotate("text", x = 1.67, y = 2.5, label = Yeast.1$Metabolites[1], size=3) +
  annotate("text", x = 1.7, y = 4, label = Yeast.1$Metabolites[2], size=3) +
  annotate("text", x = 1.7, y = 4.8, label = Yeast.1$Metabolites[3], size=3) +
  annotate("text", x = 1.7, y = 5.3, label = Yeast.1$Metabolites[4], size=3) +
  annotate("text", x = 1.7, y = 0, label = paste("Total MCM =", round(sum(Yeast$MCM, na.rm = TRUE)), "mM"), size=6) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_F)

####Make Mammal.iBMK figures####

#head(df)

#unique(df$Species)

Mammal.iBMK <- subset(df, Species == "Mammal iBMK")

#colnames(Mammal.iBMK)

#head(Mammal.iBMK)

#Calculate the Sum of the metabolites

a <- Mammal.iBMK$Concentration[1]

for (i in 2:length(Mammal.iBMK$Concentration)){
  a <- c(a, a[i-1] + Mammal.iBMK$Concentration[i])
}

Mammal.iBMK$Metabolites.sum <- a

Mammal.iBMK$Metabolites.sum

#colnames(Mammal.iBMK)

#plot(Mammal.iBMK$Concentration, Mammal.iBMK$MCM)
#abline(a = 0, b = 1)

#Calculate the Sum of the MCM

a <- Mammal.iBMK$MCM[1]

for (i in 2:length(Mammal.iBMK$Concentration)){
  if (is.na(Mammal.iBMK$MCM[i])){
    a <- c(a, a[i-1] + 0)
  }
  if (is.na(Mammal.iBMK$MCM[i]) == FALSE){
    a <- c(a, a[i-1] + Mammal.iBMK$MCM[i])
  }
}

Mammal.iBMK$MCM.sum <- a

####Identify Mg binding strength####

#range(Mammal.iBMK$Kd.app, na.rm = TRUE)

a <- c()

for (i in 1:length(Mammal.iBMK$Kd.app)){
  print(i)
  if (is.na(Mammal.iBMK$Kd.app[i])){a[i] <- "No chelation"}
  if (is.na(Mammal.iBMK$Kd.app[i]) == FALSE){
    if (Mammal.iBMK$Kd.app[i] >= 2){a[i] <- "weak"}
    if (Mammal.iBMK$Kd.app[i] < 2){a[i] <- "strong"}
  }
}

Mammal.iBMK$Mg.binding.strength <- a

####Figure G####

#colnames(Mammal.iBMK)

Figure_G <- ggplot(data = Mammal.iBMK) +
  geom_bar(mapping = aes(x = Metabolites, y = Concentration),
           stat = "identity", width = 1, fill = "#440154FF") +
  geom_bar(mapping = aes(x = Metabolites, y = MCM),
           stat = "identity", width = 1, fill = "#35B779FF") +
  scale_x_discrete(limits = Mammal.iBMK$Metabolites) +
  geom_point(aes(x= Metabolites, y= Metabolites.sum), stat="identity", color = "#440154FF")+
  geom_line(aes(x= Metabolites, y= Metabolites.sum, group = 1),
            stat="identity", color = "#440154FF", size = 1)+
  geom_point(aes(x= Metabolites, y= MCM.sum), stat="identity", color = "#35B779FF")+
  geom_line(aes(x= Metabolites, y= MCM.sum, group = 1),
            stat="identity", color = "#35B779FF", size = 1)+
  geom_hline(yintercept = sum(Mammal.iBMK$Concentration, na.rm = TRUE), color =  "#440154FF", size = 1) +
  annotate("text", x = 15, y = 250, label = paste("Estimated total metabolites = ", round(sum(Mammal.iBMK$Concentration, na.rm = TRUE))," mM", sep =""), size = 4, color = "#440154FF") +
  geom_hline(yintercept = sum(Mammal.iBMK$MCM, na.rm = TRUE), color =  "#35B779FF", size = 1) +
  annotate("text", x = 25, y = 57, label = paste("Estimated total chelated Mg2+ = ", round(sum(Mammal.iBMK$MCM, na.rm = TRUE))," mM", sep =""), size = 4, color = "#35B779FF") +
  xlim(Mammal.iBMK$Metabolites[1:40]) +
  ylim(0,250) +
  ylab("Concentration (mM)")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "Black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 12),
        legend.text = element_text(color = "Black", size = 16),
        legend.title = element_text(color = "Black", size = 16))

#ggsave("Test.png", Figure_G)


####Make Figure H####

#colnames(Mammal.iBMK)

Mammal.iBMK.1 = Mammal.iBMK %>% filter(Kd.app < 2)
Mammal.iBMK.2 = Mammal.iBMK %>% filter(Kd.app >= 2)
#View(Mammal.iBMK.2)

Figure_H <- ggplot(Mammal.iBMK, aes(x = "", y = Concentration, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  annotate("text", x = 1.75, y = 30, label = Mammal.iBMK.2$Metabolites[1], size=3) +
  annotate("text", x = 1.7, y = 70, label = Mammal.iBMK.2$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 85, label = Mammal.iBMK.2$Metabolites[3], size=3) +
  annotate("text", x = 1.8, y = 0, label = paste("Total Metabolite =", round(sum(Mammal.iBMK$Concentration)), "mM"), size=6) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_H)

####Make Figure I####

Figure_I <- ggplot(Mammal.iBMK, aes(x = "", y = MCM, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(4,6)]) +
  annotate("text", x = 1.67, y = 0.5, label = Mammal.iBMK.2$Metabolites[1], size=3) +
  annotate("text", x = 1.7, y = 1.3, label = Mammal.iBMK.2$Metabolites[4], size=3) +
  annotate("text", x = 1.67, y = 4, label = Mammal.iBMK.1$Metabolites[1], size=3) +
  annotate("text", x = 1.7, y = 7, label = Mammal.iBMK.1$Metabolites[2], size=3) +
  annotate("text", x = 1.7, y = 8, label = Mammal.iBMK.1$Metabolites[3], size=3) +
  annotate("text", x = 1.7, y = 8.75, label = Mammal.iBMK.1$Metabolites[4], size=3) +
  annotate("text", x = 1.7, y = 0, label = paste("Total MCM =", round(sum(Mammal.iBMK$MCM, na.rm = TRUE)), "mM"), size=6) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        axis.text.x=element_blank())

#ggsave("Test.png", Figure_I)

####Print figures####

#?plot_grid

Final_figure = plot_grid(plotlist = list(Figure_A, Figure_B, Figure_C,
                         Figure_D, Figure_E, Figure_F,
                         Figure_G, Figure_H, Figure_I),
                         nrow = 3)

#list.files("Read_me_figure")

ggsave("Read_me_figure/Test.svg", Final_figure, height = 24, width = 24, scale = 0.75)
