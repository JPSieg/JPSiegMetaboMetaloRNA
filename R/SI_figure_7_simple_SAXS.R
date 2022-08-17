library(tidyverse)
library(viridis)
library(cowplot)


list.files("Figures/SI_figure_7_simple_SAXS")

####Make a list of files####

vector.paths = paste("Figures/SI_figure_7_simple_SAXS/DAT_files",
      list.files("Figures/SI_figure_7_simple_SAXS/DAT_files"),
      sep = "/")

vector.P.R = paste("Figures/SI_figure_7_simple_SAXS/OUT_files",
                     list.files("Figures/SI_figure_7_simple_SAXS/OUT_files"),
                     sep = "/")

####Read in and format data####

read.dat = function(file_path = "Figures/SI_figure_7_simple_SAXS/DAT_files/25mMMg-RNA-rep1_subtr.dat"){
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

list.files("Figures/SI_figure_7_simple_SAXS/DAT_files")

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

list.files("Figures/SI_figure_7_simple_SAXS/OUT_files")

read.out = function(file_path = "Figures/SI_figure_7_simple_SAXS/OUT_files/25mMMg-RNA-rep1_subtr.out"){
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


####Fit data to WAXSiS crystal structure####

list.files("Figures/SI_figure_7_simple_SAXS/waxsresults.1474136")

Lines = readLines(file("Figures/SI_figure_7_simple_SAXS/waxsresults.1474136/intensity.dat"))

q = c()
Iq = c()
sigmaq = c()

for (i in 5:length(Lines)){
  v.line = strsplit(Lines[i], split = " ")[[1]]
  v.line = v.line[-which(v.line == "")]
  q = c(q, as.numeric(v.line[1]))
  Iq = c(Iq, as.numeric(v.line[2]))
  sigmaq = c(sigmaq, as.numeric(v.line[3]))
}

head(df)

df.model = data.frame(q, Iq, sigmaq)

plot(df.model$q, (df.model$q^2)*df.model$Iq)

list.df.fit = {}

for (i in 1:length(unique(df$Condition))){
  Condition = unique(df$Condition)[i]

  df.fit.data = df %>% filter(Condition == Condition)

  df.model.data = df.model %>%
    filter(q >= min(df.fit.data$q, na.rm = T)) %>%
    filter(q <= max(df.fit.data$q, na.rm = T))

  Iq.exp = c()
  SE.Iq = c()

  for (j in 1:nrow(df.model.data)){
    Iq.exp[j] = df.fit.data$Iq[which.min(abs(df.fit.data$q - df.model.data$q[j]))]
    SE.Iq[j] = df.fit.data$err.q[which.min(abs(df.fit.data$q - df.model.data$q[j]))]
  }

  df.model.data$Iq.exp = Iq.exp
  df.model.data$SE.Iq = SE.Iq

  fit = lm(Iq.exp ~ Iq, df.model.data)

  df.model.data$Fit.data = predict(fit)

  df.model.data$Condition = Condition

  list.df.fit[[i]] = df.model.data
}

df.fit = bind_rows(list.df.fit)

df.fit$q2.Iq = (1000)*(df.fit$q^2)*df.fit$Fit.data

####Dimensionless Kratky plot####

head(df)
head(df.fit)

Kratky.plot = ggplot() +
  #facet_wrap(~Condition) +
  geom_point(data = df,
             mapping = aes(x = q, y = q2.Iq, color = Condition),
             alpha = 0.6)  +
  geom_line(data = df.fit,
            mapping = aes(x = q, y = q2.Iq),
            color = "black") +
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

Kratky.plot

####Prod plot####

list.files("Figures/SI_figure_7_simple_SAXS")

df.out.fit = read.csv("Figures/SI_figure_7_simple_SAXS/WAXSiS_pr.csv")

Pr.plot = ggplot() +
  geom_point(data = df.out,
             mapping = aes(x = R, y = 1000*P.R, color = Condition),
             alpha = 0.5)  +
  geom_line(data = df.out.fit,
            mapping = aes(x = r, y = 0.00003*p.r.)) +
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

Pr.plot

####DENSS####

list.files("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS")

Mg.free.2mM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/2mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Mg.free.2mM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/2mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Eco80 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/Eco80_angle_0deg.png") +
  theme(panel.background = element_blank())

Eco80.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/Eco80_angle_90deg.png") +
  theme(panel.background = element_blank())

NTPCM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/NTPCM_angle_0deg.png") +
  theme(panel.background = element_blank())

NTPCM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/NTPCM_angle_90deg.png") +
  theme(panel.background = element_blank())

WMCM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/WMCM_angle_0deg.png") +
  theme(panel.background = element_blank())

WMCM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/WMCM_angle_90deg.png") +
  theme(panel.background = element_blank())

Free25 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/25mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Free25.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_WAXSiS/25mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Figure_C.DENSS = plot_grid(Mg.free.2mM, Mg.free.2mM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())


Figure_D.DENSS = plot_grid(Eco80, Eco80.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())


Figure_E.DENSS = plot_grid(NTPCM, NTPCM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

Figure_F.DENSS = plot_grid(WMCM, WMCM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

Figure_G.DENSS = plot_grid(Free25, Free25.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

####Bead####

list.files("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead")

Mg.free.2mM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/2mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Mg.free.2mM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/2mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Eco80 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/Eco80_angle_0deg.png") +
  theme(panel.background = element_blank())

Eco80.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/Eco80_angle_90deg.png") +
  theme(panel.background = element_blank())

NTPCM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/NTPCM_angle_0deg.png") +
  theme(panel.background = element_blank())

NTPCM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/NTPCM_angle_90deg.png") +
  theme(panel.background = element_blank())

WMCM = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/WMCM_angle_0deg.png") +
  theme(panel.background = element_blank())

WMCM.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/WMCM_angle_90deg.png") +
  theme(panel.background = element_blank())

Free25 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/25mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Free25.90 = ggplot() +
  draw_image("Figures/SI_figure_7_simple_SAXS/Images_from_pymol_Bead/25mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Figure_C.Bead = plot_grid(Mg.free.2mM, Mg.free.2mM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())


Figure_D.Bead = plot_grid(Eco80, Eco80.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())


Figure_E.Bead = plot_grid(NTPCM, NTPCM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

Figure_F.Bead = plot_grid(WMCM, WMCM.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

Figure_G.Bead = plot_grid(Free25, Free25.90, nrow = 1) +
  draw_image("Figures/Figure_3_degredation/Rotate_symbol.png", scale = 0.2) +
  theme(panel.background = element_blank())

####Consolidate reconstructions####

Figure_C = plot_grid(Figure_C.DENSS, ncol = 1, label_y = 0.8,
          labels = c("Electron density", "Bead model"),
          label_size = 14)
Figure_D = plot_grid(Figure_D.DENSS, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_E = plot_grid(Figure_E.DENSS, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_F = plot_grid(Figure_F.DENSS, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_G = plot_grid(Figure_G.DENSS, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)

####Consolidate plot####

SI_figure_XAB = plot_grid(Kratky.plot, Pr.plot,
                        labels = c("A", "B", "C (2 mM free Mg2+)"),
                        label_size = 20,
                        nrow = 1)

SI_Figure_XCDEFG = plot_grid(Figure_C, Figure_D, Figure_E, Figure_F, Figure_G,
                           nrow = 2,
                           labels = c("C (2 mM free)", "D (Eco80)", "E (NTPCM)", "F (WMCM)", "G (25 mM free Mg2+)"),
                           label_size = 20)

SI_figure_X = plot_grid(SI_figure_XAB,
                        SI_Figure_XCDEFG,
                        ncol = 1,
                        rel_heights = c(1, 1))

####Print plot####

ggsave("Figures/SI_figure_7_simple_SAXS/SI_figure_x.svg", SI_figure_X, scale = 1.4, height = 6, width = 7)

