library(tidyverse)
library(viridis)
library(cowplot)

####Make a list of files####

vector.paths = paste("Figures/SI_figure_x_simple_SAXS/DAT_files",
      list.files("Figures/SI_figure_x_simple_SAXS/DAT_files"),
      sep = "/")

vector.P.R = paste("Figures/SI_figure_x_simple_SAXS/OUT_files",
                     list.files("Figures/SI_figure_x_simple_SAXS/OUT_files"),
                     sep = "/")

####Read in and format data####

read.dat = function(file_path = "Figures/SI_figure_x_simple_SAXS/DAT_files/25mMMg-RNA-rep1_subtr.dat"){
  con = file(file_path)
  Lines = readLines(con)
  close(con)

  q = c()
  Iq = c()
  err.q = c()

  for (i in 9:762){
    Line.v = as.numeric(strsplit(Lines[i], split = " ")[[1]])
    q[i] = Line.v[1]
    Iq[i] = Line.v[2]
    err.q[i] = Line.v[3]
  }

  df = data.frame(q,
                  Iq,
                  err.q)

  output = df[-c(1:16),]
}

list.df = lapply(vector.paths, read.dat)

list.files("Figures/SI_figure_x_simple_SAXS/DAT_files")

list.df[[1]]$Condition = "25 mM Mg2+"
list.df[[2]]$Condition = "2 mM Mg2+"
list.df[[3]]$Condition = "NTPCM"
list.df[[4]]$Condition = "Eco80"
list.df[[5]]$Condition = "WMCM"

df = bind_rows(list.df)

df$q2.Iq = (1000)*(df$q^2)*df$Iq
df$err.q2.Iq = (1000)*(df$q^2)*df$err.q

df$Condition = factor(df$Condition,
                      levels = c("2 mM Mg2+",
                                 "NTPCM",
                                 "WMCM",
                                 "Eco80",
                                 "25 mM Mg2+"))

list.files("Figures/SI_figure_x_simple_SAXS/OUT_files")

read.out = function(file_path = "Figures/SI_figure_x_simple_SAXS/OUT_files/25mMMg-RNA-rep1_subtr.out"){
  con = file(file_path)
  Lines = readLines(con)
  close(con)

  R = c()
  P.R = c()

  for (i in 357:length(Lines)){
    Line.v = as.numeric(strsplit(Lines[i], split = " ")[[1]])
    R[i] = Line.v[3]
    P.R[i] = Line.v[5]
  }

  df = data.frame(R,
                  P.R)

  output = df[-c(1:356),]
}

list.out = lapply(vector.P.R, read.out)

list.out[[1]]$Condition = "25 mM Mg2+"
list.out[[2]]$Condition = "2 mM Mg2+"
list.out[[3]]$Condition = "NTPCM"
list.out[[4]]$Condition = "Eco80"
list.out[[5]]$Condition = "WMCM"

df.out = bind_rows(list.out)

df.out$Condition = factor(df.out$Condition,
                      levels = c("2 mM Mg2+",
                                 "NTPCM",
                                 "WMCM",
                                 "Eco80",
                                 "25 mM Mg2+"))

####Dimensionless Kratky plot####

Kratky.plot = ggplot(df, aes(x = q, y = q2.Iq, ymin = q2.Iq - err.q2.Iq, ymax = q2.Iq + err.q2.Iq, color = Condition)) +
  geom_point(alpha = 0.5)  +
  #geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "red") +
  #facet_wrap(~Condition, nrow = 1) +
  ylab(bquote(~q^2 ~I(q) (1000))) +
  xlab("q (1/\u212b)") +
  xlim(0, 0.3) +
  ylim(-2, 5) +
  theme_classic() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_text(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.25, 0.2))

####Prod plot####

Pr.plot = ggplot(df.out, aes(x = R, y = 1000*P.R, color = Condition)) +
  geom_point(alpha = 0.5)  +
  ylab("p(r) (1000)") +
  xlab("Distance (\u212b)") +
  theme_classic() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_text(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.35, 0.3))

####Fit comparisons####

list.files("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits")

read.waxsis = function(data_path){

  con = file(data_path)
  Lines = readLines(con)
  close(con)

  q = c()
  Iq.fit = c()
  Iq.fit.err = c()
  Iq.calc = c()
  Iq.calc.err = c()

  for (i in 10:length(Lines)){
    Line.vector = strsplit(Lines[[i]], split = " ")[[1]]
    Line.vector = Line.vector[-which(Line.vector == "")]
    Line.vector = as.numeric(Line.vector)
    q[i] = Line.vector[1]
    Iq.fit[i] = Line.vector[2]
    Iq.fit.err[i] = Line.vector[3]
    Iq.calc[i] = Line.vector[4]
    Iq.calc.err[i] = Line.vector[5]
  }

  output = data.frame(q,
                         Iq.fit,
                         Iq.fit.err,
                         Iq.calc,
                         Iq.calc.err)

}

df.2mM = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/2mMMg_over.xvg")
df.2mM$Condition = "2 mM Mg X2 = 1.09"
df.2mM$State = "Monomer"
df.25mM = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/25mMMg_over.xvg")
df.25mM$Condition = "25 mM Mg X2 = 2.7"
df.25mM$State = "Monomer"
df.WMCM = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/WMCM_over.xvg")
df.WMCM$Condition = "WMCM X2 = 1.534"
df.WMCM$State = "Monomer"
df.25mM$State = "Monomer"
df.NTPCM = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/NTPCM_over.xvg")
df.NTPCM$Condition = "NTPCM X2 = 0.86"
df.NTPCM$State = "Monomer"
df.Eco80 = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/Eco80_over.xvg")
df.Eco80$Condition = "Eco80 X2 = 0.81"
df.Eco80$State = "Monomer"

df.2mM.dimer = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/2mM_dimer_over.xvg")
df.2mM.dimer$Condition = "2 mM Mg X2 = 1.09"
df.2mM.dimer$State = "Dimer"
df.25mM.dimer = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/25mM_dimer_over.xvg")
df.25mM.dimer$Condition = "25 mM Mg X2 = 2.7"
df.25mM.dimer$State = "Dimer"
df.WMCM.dimer = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/WMCM_dimer_over.xvg")
df.WMCM.dimer$Condition = "WMCM X2 = 1.534"
df.WMCM.dimer$State = "Dimer"
df.25mM.dimer$State = "Dimer"
df.NTPCM.dimer = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/NTPCM_dimer_over.xvg")
df.NTPCM.dimer$Condition = "NTPCM X2 = 0.86"
df.NTPCM.dimer$State = "Dimer"
df.Eco80.dimer = read.waxsis("Figures/SI_figure_x_simple_SAXS/WAXSiS_fits/Eco80_dimer_over.xvg")
df.Eco80.dimer$Condition = "Eco80 X2 = 0.81"
df.Eco80.dimer$State = "Dimer"

df.exp = bind_rows(df.2mM, df.25mM, df.NTPCM, df.WMCM, df.Eco80,
                   df.2mM.dimer, df.25mM.dimer, df.NTPCM.dimer, df.WMCM.dimer, df.Eco80.dimer)

df.exp$Condition = factor(df.exp$Condition,
                          levels = c("2 mM Mg X2 = 1.09",
                                     "NTPCM X2 = 0.86",
                                     "WMCM X2 = 1.534",
                                    "Eco80 X2 = 0.81",
                                    "25 mM Mg X2 = 2.7"))

head(df.exp)
unique(df.exp$State)
df.exp[which(is.na(df.exp$State)),]

Residuals = ggplot(df.exp, aes(x = q, y = (Iq.fit - Iq.calc)/Iq.fit, color = Condition, group = State)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(limits = c(-5, 5), n.breaks = 3) +
  geom_hline(yintercept = 0) +
  ylab("Normalized\nResiduals") +
  theme_classic() +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  facet_wrap(c("Condition", "State")) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 12),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_text(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = "none")

Data = ggplot() +
  geom_point(data = df.exp, mapping = aes(x = q, y = Iq.fit, color = Condition), alpha = 0.5) +
  geom_line(data = df.exp, mapping = aes(x = q, y = Iq.calc, shape = State, group = State)) +
  scale_color_manual(values = c("dimgrey", viridis(n =  7)[c(3, 1, 6)], "red")) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~Condition, nrow = 1) +
  xlab("q (1/\u212b)") +
  ylab("I(q)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 16),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        strip.background = element_rect(size = 1),
        strip.text = element_text(color = "Black", size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = "none")

Data

SI_figure_XC = plot_grid(Residuals, Data, ncol = 1, rel_heights = c(1.5, 4),
          align = "v")

####Consolidate plot####

SI_figure_XAB = plot_grid(Kratky.plot, Pr.plot,
                        labels = c("A", "B"),
                        label_size = 20)

SI_figure_X = plot_grid(SI_figure_XAB,
                        SI_figure_XC,
                        ncol = 1,
                        labels = c("", "C"),
                        label_size = 20,
                        rel_heights = c(2, 3))

####Print plot####

ggsave("Figures/SI_figure_x_simple_SAXS/SI_figure_x.png", SI_figure_X, scale = 1.4, height = 7, width = 15)

