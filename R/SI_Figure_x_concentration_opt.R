library(MeltR)
library(tidyverse)
library(viridis)
library(cowplot)

####Make error combinations 100*100 = 10,000#####

error = seq(-50, 50, length.out = 100)

list.df = {}

for (i in 1:100){
  list.df[[i]] = data.frame("FAM.error" = error[i],
                            "BHQ1.error" = error)
}

df.model = bind_rows(list.df)

#(df.model)

####Model data no concentration error####

print("Fitting model data no concentration error")

#list.files("Figures")

#df.no.error = meltR.F.model(-56.2, -0.1364, Emission_SD = 0.05)

#write.csv(df.no.error, "Figures/SI_figure_x_concentration_correction/Model_data_no_conc_error.csv", row.names = FALSE)


df.error = meltR.F.model(-56.2, -0.1364, FAM_error = 1.2, BHQ1_error = 0.8, Emission_SD = 0.05)

pb = txtProgressBar(min = 1, max = nrow(df.model), initial = 1)

dG = c()
dH = c()
dS = c()

for (i in 1:nrow(df.model)){
  setTxtProgressBar(pb, i)
  df = read.csv("Figures/SI_figure_x_concentration_correction/Model_data_no_conc_error.csv")
  tryCatch({df$A = df$A*(df.model$FAM.error[i]/100 + 1)
  df$B = df$B*(df.model$BHQ1.error[i]/100 + 1)
  fit = meltR.F(df, Optimize_conc = FALSE,
                K_range = c(10, 200),
                K_error_quantile = 1,
                silent = TRUE)
  dG[i] = fit$VantHoff$G[1]
  dH[i] = fit$VantHoff$H[1]
  dS[i] = fit$VantHoff$S[1]},
           error = function(e){dG[i] = NA
           dH[i] = NA
           dS[i] = NA})


  close(pb)
}

dG.known = -56.2 + (273.15+37)*0.1364
df.model$dH.error = abs(dH - -56.2)
df.model$dS.error = abs(dS - -136.4)
df.model$dG.error = abs(dG - dG.known)



####Model data with concentration error####

print("Fitting model data with concentration error ")

#list.files("Figures")

#df.error = meltR.F.model(-56.2, -0.1364, FAM_error = 1.15, Emission_SD = 0.05)

#write.csv(df.error, "Figures/SI_figure_x_concentration_correction/Model_data_with_conc_error.csv", row.names = FALSE)

pb = txtProgressBar(min = 1, max = nrow(df.model), initial = 1)

dG = c()
dH = c()
dS = c()

for (i in 1:nrow(df.model)){
  setTxtProgressBar(pb, i)
  df = read.csv("Figures/SI_figure_x_concentration_correction/Model_data_with_conc_error.csv")
  tryCatch({df$A = 200*(df.model$FAM.error[i]/100 + 1)
  df$B = df$B*(df.model$BHQ1.error[i]/100 + 1)
  fit = meltR.F(df, Optimize_conc = FALSE,
                K_range = c(10, 200),
                K_error_quantile = 1,
                silent = TRUE)
  dG[i] = fit$VantHoff$G[1]
  dH[i] = fit$VantHoff$H[1]
  dS[i] = fit$VantHoff$S[1]},
  error = function(e){dG[i] = NA
  dH[i] = NA
  dS[i] = NA})


  close(pb)
}

df.model.error = df.model
df.model.error$dH.error = abs(dH - -56.2)
df.model.error$dS.error = abs(dS - -136.4)
df.model.error$dG.error = abs(dG - dG.known)

####Model data with concentration error and optimization on####

print("Fitting model data with concentration error and optimization on")

dG = c()
dH = c()
dS = c()

for (i in 1:nrow(df.model)){
  setTxtProgressBar(pb, i)
  df = read.csv("Figures/SI_figure_x_concentration_correction/Model_data_with_conc_error.csv")
  tryCatch({df$A = 200*(df.model$FAM.error[i]/100 + 1)
  df$B = df$B*(df.model$BHQ1.error[i]/100 + 1)
  fit = meltR.F(df, Optimize_conc = TRUE,
                K_range = c(10, 200),
                K_error_quantile = 1,
                silent = TRUE)
  dG[i] = fit$VantHoff$G[1]
  dH[i] = fit$VantHoff$H[1]
  dS[i] = fit$VantHoff$S[1]},
  error = function(e){dG[i] = NA
  dH[i] = NA
  dS[i] = NA})


  close(pb)
}

df.model.optimize = df.model
df.model.optimize$dH.error = abs(dH - -56.2)
df.model.optimize$dS.error = abs(dS - -136.4)
df.model.optimize$dG.error = abs(dG - dG.known)


####Make plots####

df.model$Error.type = "None"
df.model.error$Error.type = "+15% FAM error"
df.model.optimize$Error.type = "+15% FAM error & optimized"
df.final = bind_rows(df.model, df.model.error, df.model.optimize)

df.final$R = (4 + 4*df.final$FAM.error/100)/(4 + 4*df.final$BHQ1.error/100)
df.final$Total = (200 + 200*df.final$FAM.error/100) + (200 + 200*df.final$BHQ1.error/100)


write.csv(df.final, "Figures/SI_figure_x_concentration_correction/Final_conc_error_10000_fits.csv", row.names = FALSE)

df.final = read.csv("Figures/SI_figure_x_concentration_correction/Final_conc_error_10000_fits.csv")

list.files()

df.final$Error.type = factor(df.final$Error.type,
                             levels = c("None", "+15% FAM error", "+15% FAM error & optimized"),
                             labels = c("None", "+20% FAM error", "+20% FAM error & optimized"))

dG.plot = ggplot(df.final, aes(x = FAM.error, y = BHQ1.error, fill = dG.error)) +
  geom_tile() +
  scale_fill_viridis(option = "F", direction = -1) +
  geom_abline(slope = 1, intercept = 1) +
  theme_classic() +
  xlab("%FAM error") +
  ylab("%BHQ1 error") +
  facet_wrap(~Error.type, nrow = 1) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))

R.plot = ggplot(df.final, aes(x = R, y = dG.error)) +
  geom_smooth() +
  geom_point() +
  scale_fill_viridis(option = "F", direction = -1) +
  geom_vline(xintercept = c(1, (4+4*0.15)/4)) +
  theme_classic() +
  scale_x_continuous(trans = "log2") +
  xlab("R = [FAM-RNA]/[RNA-BHQ1]") +
  ylab("dG error (kcal/mol)") +
  facet_wrap(~Error.type, ncol = 1) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))

Total.plot = ggplot(df.final, aes(x = Total, y = dG.error)) +
  geom_point() +
  scale_fill_viridis(option = "F", direction = -1) +
  theme_classic() +
  xlab("Total concentration error (nm)") +
  ylab("dG error (kcal/mol)") +
  facet_wrap(~Error.type, ncol = 1) +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))

dH.plot = ggplot(df.final, aes(x = FAM.error, y = BHQ1.error, fill = dH.error)) +
  geom_tile() +
  facet_wrap(~Error.type, ncol = 1) +
  scale_fill_viridis(option = "G", direction = -1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  xlab("%FAM error") +
  ylab("%BHQ1 error") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))


dS.plot = ggplot(df.final, aes(x = FAM.error, y = BHQ1.error, fill = dS.error)) +
  geom_tile() +
  facet_wrap(~Error.type, ncol = 1) +
  scale_fill_viridis(option = "F", direction = -1) +
  geom_abline(slope = 1, intercept = 1) +
  theme_classic() +
  xlab("%FAM error") +
  ylab("%BHQ1 error") +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))

list.files("Figures/SI_figure_x_concentration_correction")


df = read.csv("Figures/SI_figure_x_concentration_correction/Model_data_with_conc_error.csv")

head(df)

Mmodel <- function(x){ (2-0.2)*(1 - ((0.1+200+x)-(((0.1+200+x)^2)-(4*200*x))^(1/2))/(2*200)) + 0.2 }
Mmodel.cor <- function(x){ (2-0.2*1.2)*(1 - ((0.1+200*1.2+x)-(((0.1+200*1.2+x)^2)-(4*200*1.2*x))^(1/2))/(2*200*1.2)) + 0.2 }

opt.plot = ggplot(df %>% filter(Reading == 1), aes(x = B, Emission)) +
  geom_point() +
  geom_function(fun = Mmodel) +
  geom_function(fun = Mmodel.cor, color = "red") +
  xlab("[BHQ1] (nM)") +
  ylab("Emission") +
  annotate("text", x = 500, y = 0.5, label = "Start R = 1") +
  annotate("text", x = 500, y = 0.4, label = "optimized R = 1.2", color = "red") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "Black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "Black", size = 16),
        axis.title.x = element_text(color = "Black", size = 16),
        axis.title.y = element_text(color = "Black", size = 16),
        legend.text = element_text(color = "Black", size = 10),
        legend.title = element_text(color = "Black", size = 10))




plot.final = plot_grid(dG.plot,
                       opt.plot,
          nrow = 1,
          labels = c("A", "B"),
          rel_widths = c(3,1))

ggsave("Figures/SI_figure_x_concentration_correction/SI_figure_x_concentration.png", width = 7,
       height = 2,
       scale = 2,
       plot.final)

df.er = df.final %>%
  filter(!is.na(dG.error)) %>%
  filter(Error.type == "+20% FAM error & optimized") %>%
  select(dG.error)

str(df.er$dG.error)

range(df.er$dG.error)
mean(df.er$dG.error)
hist(df.er$dG.error)
quantile(df.er$dG.error, 0.80)
