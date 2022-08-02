library(MeltR)
library(tidyverse)
library(viridis)
library(cowplot)








df = meltR.F.model(-56.2, -0.1364, Emission_SD = 0)

df$Kd = (10^9)/exp((df$S/0.0019872)-(df$H/((273.15+df$Temperature)*0.0019872)))

df$Kd.to.A = df$Kd/df$A

range(df$Kd.to.A)
range(df$Kd)
df$Kd.range = NA

df$Kd.range[which(df$Kd.to.A <= 0.1)] = "KD < [FAM-RNA]/10"
df$Kd.range[which(df$Kd.to.A >= 10)] = "KD > 10x[FAM-RNA]"
df$Kd.range[which(is.na(df$Kd.range))] = "Optimum range"

range(df$Kd.to.A)
head(df)

ggplot(df , aes(x = B, y = Emission, group = Reading, color = log10(Kd))) +
  geom_line() +
  facet_wrap(~Kd.range) +
  scale_color_viridis(option = "H", name = "log10(KD)") +
  theme_classic() +
  xlab("[RNA-BHQ1] (nM)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.text = element_markdown(color = "Black", size = 14),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = c(0.5, 0.35))

ggsave("Figures/SI_Figure_x_optimal_Kd_range/Optimal_Kd_range.png", scale = 2)


