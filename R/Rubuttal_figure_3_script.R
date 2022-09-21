library(tidyverse)

ddG.calc = function(free.Mg = 0.001, N = 8, Na = (0.240 + 0.140)){

  ####At Mg concentration i####

  X1.i = Na/(Na + ((8.1 - 32.4/N)*(5.2-log(Na))*free.Mg))

  a1.i = -0.075*log(Na) + 0.012*(log(Na)^2)
  b1.i = 0.018*(log(Na)^2)
  dg1.i = a1.i + b1.i/N

  a2.i = (-6/N) + 0.025*log(free.Mg) + 0.0068*(log(free.Mg)^2)
  b2.i = log(free.Mg) + 0.38*(log(free.Mg)^2)
  dg2.i = a2.i + b2.i/(N^2)

  X2.i = 1 - X1.i

  dg12.i = -0.6*X1.i*X2.i*log(Na)*log(((1/X1.i)-1)*Na)/N

  dG.i = (N-1)*(X1.i*dg1.i + X2.i*dg2.i) + dg12.i

  ####At Mg concentration 2####

  X1.2 = Na/(Na + ((8.1 - 32.4/N)*(5.2-log(Na))*0.002))

  a1.2 = -0.075*log(Na) + 0.012*(log(Na)^2)
  b1.2 = 0.018*(log(Na)^2)
  dg1.2 = a1.2 + b1.2/N

  a2.2 = (-6/N) + 0.025*log(0.002) + 0.0068*(log(0.002)^2)
  b2.2 = log(0.002) + 0.38*(log(0.002)^2)
  dg2.2 = a2.2 + b2.2/(N^2)

  X2.2 = 1 - X1.2

  dg12.2 = -0.6*X1.2*X2.2*log(Na)*log(((1/X1.2)-1)*Na)/N

  dG.2 = (N-1)*(X1.2*dg1.2 + X2.2*dg2.2) + dg12.2

  ddG = dG.i - dG.2
}

log10.Mg.free = seq(log10(0.0001), log10(0.1), length.out = 500)
Mg.free = 10^log10.Mg.free
ddG = ddG.calc(Mg.free)

df = data.frame(Mg.free, ddG)

ggplot() +
  geom_polygon(mapping = aes(x = c(1.3, 3.5, 3.5, 1.3), y = c(-0.6, -0.6, 0.6, 0.6)),
               alpha = 0.25, fill = "blue") +
  geom_polygon(mapping = aes(x = c(1.8, 2.2, 2.2, 1.8), y = c(-0.6, -0.6, 0.6, 0.6)),
               alpha = 0.5, fill = "blue") +
  geom_hline(yintercept = c(-0.25, 0, 0.25)) +
  geom_vline(xintercept = 2, color = "blue") +
  geom_point(data = df, mapping = aes(x = 1000*Mg.free, y = ddG)) +
  scale_x_continuous(limits = c(0.1, 4)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, -0.25, 0, 0.25, 0.5), expand = c(0, 0)) +
  theme_classic() +
  xlab(bquote('Free'~'Mg'^'2+'~'(mM)')) +
  ylab("\u0394\u0394G\u00B037 (kcal/mol)") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.76),
        legend.background = element_blank(),
        plot.title = element_text(color = "Black", size = 14,hjust = 0.5))

list.files("Reviewer_response_figures")

ggsave("Reviewer_response_figures/Rebutal_figure_3.svg", scale = 1.5)
