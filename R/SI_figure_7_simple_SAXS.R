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


####DENSS####

list.files("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS")

Mg.free.2mM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/2mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Mg.free.2mM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/2mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Eco80 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/Eco80_angle_0deg.png") +
  theme(panel.background = element_blank())

Eco80.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/Eco80_angle_90deg.png") +
  theme(panel.background = element_blank())

NTPCM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/NTPCM_angle_0deg.png") +
  theme(panel.background = element_blank())

NTPCM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/NTPCM_angle_90deg.png") +
  theme(panel.background = element_blank())

WMCM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/WMCM_angle_0deg.png") +
  theme(panel.background = element_blank())

WMCM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/WMCM_angle_90deg.png") +
  theme(panel.background = element_blank())

Free25 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/25mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Free25.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_WAXSiS/25mMMg_angle_90deg.png") +
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

list.files("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead")

Mg.free.2mM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/2mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Mg.free.2mM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/2mMMg_angle_90deg.png") +
  theme(panel.background = element_blank())

Eco80 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/Eco80_angle_0deg.png") +
  theme(panel.background = element_blank())

Eco80.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/Eco80_angle_90deg.png") +
  theme(panel.background = element_blank())

NTPCM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/NTPCM_angle_0deg.png") +
  theme(panel.background = element_blank())

NTPCM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/NTPCM_angle_90deg.png") +
  theme(panel.background = element_blank())

WMCM = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/WMCM_angle_0deg.png") +
  theme(panel.background = element_blank())

WMCM.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/WMCM_angle_90deg.png") +
  theme(panel.background = element_blank())

Free25 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/25mMMg_angle_0deg.png") +
  theme(panel.background = element_blank())

Free25.90 = ggplot() +
  draw_image("Figures/SI_figure_x_simple_SAXS/Images_from_pymol_Bead/25mMMg_angle_90deg.png") +
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

Figure_C = plot_grid(Figure_C.DENSS, Figure_C.Bead, ncol = 1, label_y = 0.8,
          labels = c("Electron density", "Bead model"),
          label_size = 14)
Figure_D = plot_grid(Figure_D.DENSS, Figure_D.Bead, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_E = plot_grid(Figure_E.DENSS, Figure_E.Bead, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_F = plot_grid(Figure_F.DENSS, Figure_F.Bead, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)
Figure_G = plot_grid(Figure_G.DENSS, Figure_G.Bead, ncol = 1, label_y = 0.8,
                     labels = c("Electron density", "Bead model"),
                     label_size = 14)

####Consolidate plot####

SI_figure_XABC = plot_grid(Kratky.plot, Pr.plot, Figure_C,
                        labels = c("A", "B", "C (2 mM free Mg2+)"),
                        label_size = 20,
                        nrow = 1,
                        rel_widths = c(24, 24, 16))

SI_Figure_XDEFG = plot_grid(Figure_D, Figure_E, Figure_F, Figure_G,
                           nrow = 1,
                           labels = c("D (Eco80)", "E (NTPCM)", "F (WMCM)", "G (25 mM free Mg2+)"),
                           label_size = 20)

SI_figure_X = plot_grid(SI_figure_XABC,
                        SI_Figure_XDEFG,
                        ncol = 1,
                        rel_heights = c(1, 1))

####Print plot####

ggsave("Figures/SI_figure_x_simple_SAXS/SI_figure_x.svg", SI_figure_X, scale = 1.4, height = 6, width = 10)

