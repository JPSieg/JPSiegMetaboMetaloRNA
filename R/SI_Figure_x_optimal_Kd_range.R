library(MeltR)
library(tidyverse)
library(viridis)
library(cowplot)

df = read.csv("Figures/SI_figure_x_concentration_correction/Model_data_no_conc_error.csv")

fit = meltR.F(df, Optimize_conc = FALSE,
                     K_range = c(0.01, 100000),
                     K_error_quantile = 1,
                     silent = TRUE)

fit$K$Kd = 10^9/fit$K$K
fit$K$lnKd = log(fit$K$Kd)
fit$K$SE.lnK = fit$K$SE.K/fit$K$K

length(unique(df$Temperature))

unique(df$Reading)


fit$K$Kd.range = NA
fit$K$Kd.range[which(fit$K$Kd <= 5)] = "Too low"
fit$K$Kd.range[which(fit$K$Kd >= 500)] = "Too high"
fit$K$Kd.range[which(is.na(fit$K$Kd.range))] = "Accurate"


fit$K$Kd.range = factor(fit$K$Kd.range,
                        levels = c("Too low", "Accurate", "Too high"))

list.df = {}

for (i in 1:nrow(fit$K)){
  df.reading = df %>% filter(Reading == (1:nrow(fit$K))[i])
  df.reading$Kd.range = fit$K$Kd.range[(1:nrow(fit$K))[i]]
  list.df[[i]] = df.reading
}

df = bind_rows(list.df)

df$Kd.range = factor(df$Kd.range,
                        levels = c("Too low", "Accurate", "Too high"))


min(df %>% filter(Kd.range == "Accurate") %>% select(Temperature))

Plot.B = ggplot(df %>% filter(Temperature != 23), aes(x = B, y = Emission, color = Temperature))+
  geom_point() +
  facet_wrap(~Kd.range, nrow = 1) +
  scale_color_viridis(option = "H") +
  theme_classic() +
  xlab("[RNA-BHQ1] (nM)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10),
        legend.position = c(0.85, 0.18))
Plot.B


Plot.A = ggplot(fit$K, aes(x = Kd, y = SE.lnK, color = Kd.range)) +
  geom_point() +
  scale_x_continuous(trans = "log10", limits = c(0.001, 10000)) +
  ylim(0, 10) +
  xlab("Kd (nM)") +
  ylab("Standard error in ln[Kd]") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.75))

SI_fig = plot_grid(Plot.A, Plot.B, nrow = 1, rel_widths = c(1,2.5), labels = c("A", "B"))

ggsave("Figures/SI_Figure_x_optimal_Kd_range/Optimal_Kd_range.png", SI_fig, width = 15)


