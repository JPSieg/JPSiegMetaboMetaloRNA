
library(tidyverse)
library(ggpubr)
library(ggtext)

list.files("Figures/SI_figure_4_AU_destabilization_versus_AU_content")

df = read.csv("Figures/SI_figure_4_AU_destabilization_versus_AU_content/ddG_data.csv")

head(df)

View(df)

ggplot(df %>% filter(Condition != "2 mM free"), aes(x = AU.content, y = ddG.kcal.mol,
                                                ymax = ddG.kcal.mol + error, ymin = ddG.kcal.mol - error,
                                                label = Label)) +
  facet_wrap(~Condition) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..), label.y = 0.7) +
  theme_classic() +
  scale_x_continuous(limits = c(20, 80), breaks = c(25, 50, 75)) +
  xlab("AU content (%)") +
  ylab("ddG (kcal/mol)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_markdown(color = "Black", size = 14),
        legend.title = element_text(color = "Black", size = 14))

ggsave("Figures/SI_figure_4_AU_destabilization_versus_AU_content/SI_figure_x_AU_destabilization_versus_AU_content.svg",
      width = 4,
       height = 2,
       scale = 1.5)
