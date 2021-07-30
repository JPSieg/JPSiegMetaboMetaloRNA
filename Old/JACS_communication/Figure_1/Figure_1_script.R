setwd("~/Jacob/Research/Manuscripts/JACS_communication/Figure_1")

####Load dependent packages####

library(tidyverse)
library(cowplot)
library(viridis)

####Load in data####

#list.files("Read_me_figure")

E.coli <- read.csv("E.coli_metabolites.csv")

####Make E.coli figures####

#head(df)

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


####Identify Mg binding strength####

#range(E.coli$Kd.app, na.rm = TRUE)

a <- c()

for (i in 1:length(E.coli$Kd)){
  print(i)
  if (is.na(E.coli$Kd[i])){a[i] <- "No chelation"}
  if (is.na(E.coli$Kd[i]) == FALSE){
    if (E.coli$Kd[i] >= 2){a[i] <- "weak"}
    if (E.coli$Kd[i] < 2){a[i] <- "strong"}
  }
}

E.coli$Mg.binding.strength <- a

####Figure A####

#colnames(E.coli)

Figure_A <- ggplot(data = E.coli) +
  geom_bar(mapping = aes(x = Metabolites, y = Concentration, fill = Mg.binding.strength),
           stat = "identity", width = 1) +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  geom_point(aes(x= Metabolites, y= Metabolites.sum), stat="identity", color = "black")+
  geom_line(aes(x= Metabolites, y= Metabolites.sum, group = 1),
            stat="identity", color = "black", size = 1)+
  geom_hline(yintercept = sum(E.coli$Concentration, na.rm = TRUE), color =  "black", size = 1) +
  annotate("text", x = 10, y = 250, label = paste("Total metabolites = ", round(sum(E.coli$Concentration, na.rm = TRUE))," mM", sep =""), size = 4, color = "Black") +
  geom_hline(yintercept = 0.8*sum(E.coli$Concentration, na.rm = TRUE), color =  "black", size = 1) +
  annotate("text", x = 10, y = 210, label = paste("85% of total metabolites = ", round(0.8*sum(E.coli$Concentration, na.rm = TRUE))," mM", sep =""), size = 4, color = "Black") +
  xlim(E.coli$Metabolites[1:length(which(E.coli$Metabolites.sum <= 0.8*sum(E.coli$Concentration)))]) +
  ylim(0,250) +
  ylab("Concentration (mM)")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text.x = element_text(color = "Black", size = 8,
                                   angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.8, 0.5),
        axis.text.y = element_text(color = "Black", size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 12),
        legend.text = element_text(color = "Black", size = 16))

#ggsave("Test.png", Figure_C)

####Write the top 14 80% of metabolites####

E.coli.3 = E.coli %>% filter(Metabolites.sum <= 0.8*sum(E.coli$Concentration))

colnames(E.coli.3)

write.csv(E.coli.3 %>% select(Metabolites, Concentration, Kd, Metabolites.sum, Mg.binding.strength),
          "Top_14_E.coli_metabolites.csv", row.names = FALSE)

####Add in ITC Kds at 37 C and read back in for circle plots####

E.coli.3 = read.csv("Top_14_E.coli_metabolites_edited.csv")


####Make Figure B####

#colnames(E.coli)

E.coli.1 = E.coli.3 %>% filter(Kd < 2)
E.coli.2 = E.coli.3 %>% filter(Kd >= 2)
#View(E.coli.2)

Figure_B <- ggplot(E.coli.3, aes(x = "", y = Concentration, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  annotate("text", x = 1.75, y = 30, label = E.coli.2$Metabolites[1], size=3) +
  annotate("text", x = 1.63, y = 95, label = E.coli.2$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 125, label = E.coli.2$Metabolites[3], size=3) +
  annotate("text", x = 1.65, y = 130, label = E.coli.2$Metabolites[4], size=3) +
  annotate("text", x = 1.65, y = 140, label = E.coli.2$Metabolites[5], size=3) +
  annotate("text", x = 1.65, y = 160, label = E.coli.1$Metabolites[1], size=3) +
  annotate("text", x = 1.65, y = 165, label = E.coli.1$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 170, label = E.coli.1$Metabolites[3], size=3) +
  annotate("text", x = 1.65, y = 175, label = E.coli.1$Metabolites[4], size=3) +
  annotate("text", x = 1.65, y = 180, label = E.coli.1$Metabolites[5], size=3) +
  annotate("text", x = 1, y = 25, label = paste(round(sum(E.coli.2$Concentration)), "mM"), size=4, color = "white") +
  annotate("text", x = 1, y = 175, label = paste(round(sum(E.coli.1$Concentration)), "mM"), size=4, color = "white") +
  annotate("text", x = 1.8, y = 0, label = paste("Metabolites =", round(0.8*sum(E.coli$Concentration)), "mM"), size=4) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.text.x=element_blank())

####Make Figure C####

E.coli.1$MCM = E.coli.1$Concentration*2/(E.coli.1$Kd + 2)
E.coli.2$MCM = E.coli.2$Concentration*2/(E.coli.2$Kd + 2)
E.coli.3$MCM = E.coli.3$Concentration*2/(E.coli.3$Kd + 2)

Figure_C <- ggplot(E.coli.3, aes(x = "", y = MCM, fill = Mg.binding.strength)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = viridis(n =  7, direction = -1)[c(2,4,6)]) +
  annotate("text", x = 1.75, y = 1.5, label = E.coli.2$Metabolites[1], size=3) +
  annotate("text", x = 1.63, y = 3.5, label = E.coli.2$Metabolites[4], size=3) +
  annotate("text", x = 1.65, y = 5, label = E.coli.2$Metabolites[5], size=3) +
  annotate("text", x = 1.65, y = 7, label = E.coli.1$Metabolites[1], size=3) +
  annotate("text", x = 1.65, y = 10, label = E.coli.1$Metabolites[2], size=3) +
  annotate("text", x = 1.65, y = 13, label = E.coli.1$Metabolites[3], size=3) +
  annotate("text", x = 1.65, y = 16, label = E.coli.1$Metabolites[4], size=3) +
  annotate("text", x = 1.65, y = 20, label = E.coli.1$Metabolites[5], size=3) +
  annotate("text", x = 1, y = 3, label = paste(round(sum(E.coli.2$MCM)), "mM"), size=4, color = "white") +
  annotate("text", x = 1, y = 16, label = paste(round(sum(E.coli.1$MCM)), "mM"), size=4, color = "white") +
  annotate("text", x = 1.8, y = 0, label = paste("MCM =", round(sum(E.coli.3$MCM)), "mM"), size=4) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(color="Black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.text.x=element_blank())


####Print figures####

?plot_grid

Figure_B_C = plot_grid(plotlist = list(Figure_B, Figure_C), labels = c("B", "C"))

Final_figure = plot_grid(plotlist = list(Figure_A, Figure_B_C),
                         ncol = 1, labels = c("A"))

#list.files("Read_me_figure")

?ggsave
ggsave("Test.svg", Final_figure, width = 3.3, height = 3.3, units = "in", scale = 2)


?cowplot
