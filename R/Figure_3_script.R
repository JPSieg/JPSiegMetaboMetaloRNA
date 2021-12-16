library(tidyverse)
library(MeltR)
library(viridis)
library(ggbeeswarm)

####Format data####

correction = function(x){ #function to apply a temp correction
  x - 0.05424*x + 4.44724
}

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")
paste("Figures/Figure_3/Index_files", list.files("Figures/Figure_3/Index_files"), sep = "/")

?AB.qPCR.export.to.MeltR.csv

#js5057 NTP-Mg F

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_F_and_G_data.txt",
                            "Figures/Figure_3/Index_files/js5057_Helix_F_index.csv",
                            "Figures/Figure_3/Fluorescence_data/js5057_Helix_F_formatted.csv",
                            correction,
                            Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5057 NTP-Mg G

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_F_and_G_data.txt",
                                 "Figures/Figure_3/Index_files/js5057_Helix_G_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5057_Helix_G_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5057 NTP-Mg H

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_H_and_I_data.txt",
                                 "Figures/Figure_3/Index_files/js5057_Helix_H_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5057_Helix_H_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5057 NTP-Mg I

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_H_and_I_data.txt",
                                 "Figures/Figure_3/Index_files/js5057_Helix_I_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5057_Helix_I_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5057 NTP-Mg J

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5057_js5060_Helix_J_data.txt",
                                 "Figures/Figure_3/Index_files/js5057_Helix_J_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5057_Helix_J_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5060 2mMMg F

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_F_and_G_data.txt",
                                 "Figures/Figure_3/Index_files/js5060_Helix_F_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5060_Helix_F_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5060 2mMMg G

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_F_and_G_data.txt",
                                 "Figures/Figure_3/Index_files/js5060_Helix_G_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5060_Helix_G_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5060 2mMMg H

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_H_and_I_data.txt",
                                 "Figures/Figure_3/Index_files/js5060_Helix_H_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5060_Helix_H_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5060 2mMMg I

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_H_and_I_data.txt",
                                 "Figures/Figure_3/Index_files/js5060_Helix_I_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5060_Helix_I_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

#js5060 2mMMg J

df = AB.qPCR.export.to.MeltR.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_J_data.txt",
                                 "Figures/Figure_3/Index_files/js5060_Helix_J_index.csv",
                                 "Figures/Figure_3/Fluorescence_data/js5060_Helix_J_formatted.csv",
                                 correction,
                                 Scale = 100000)

df = df %>% filter(Reading == 1)

plot(df$B, df$Emission) #looks good

####Fit data####

list.df = {}
list.fit = {}

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

#NTPCM Helix F

list.df[[1]] = read.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_F_formatted.csv")

?meltR.F

list.fit[[1]] = meltR.F(list.df[[1]],
        file_prefix = "Figures/Figure_3/Fluorescence_data/16mMNTPCM_Helix_F",
        Save_results = "all",
        B.conc.Tm = 175)

list.fit[[1]]$Tms$Helix = "F"
list.fit[[1]]$Tms$Condition = "16mMMg-NTPs"
list.fit[[1]]$VantHoff$Helix = "F"
list.fit[[1]]$VantHoff$Condition = "16mMMg-NTPs"

#NTPCM Helix G

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[2]] = read.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_G_formatted.csv")

list.fit[[2]] = meltR.F(list.df[[2]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/16mMNTPCM_Helix_G",
                        Save_results = "all",
                        B.conc.Tm = 250)

list.fit[[2]]$Tms$Helix = "G"
list.fit[[2]]$Tms$Condition = "16mMMg-NTPs"
list.fit[[2]]$VantHoff$Helix = "G"
list.fit[[2]]$VantHoff$Condition = "16mMMg-NTPs"

#NTPCM Helix H

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[3]] = read.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_H_formatted.csv")

list.fit[[3]] = meltR.F(list.df[[3]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/16mMNTPCM_Helix_H",
                        Save_results = "all",
                        B.conc.Tm = 300)

list.fit[[3]]$Tms$Helix = "H"
list.fit[[3]]$Tms$Condition = "16mMMg-NTPs"
list.fit[[3]]$VantHoff$Helix = "H"
list.fit[[3]]$VantHoff$Condition = "16mMMg-NTPs"

#NTPCM Helix I

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[4]] = read.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_I_formatted.csv")

list.fit[[4]] = meltR.F(list.df[[4]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/16mMNTPCM_Helix_I",
                        Save_results = "all",
                        B.conc.Tm = 200)

list.fit[[4]]$Tms$Helix = "I"
list.fit[[4]]$Tms$Condition = "16mMMg-NTPs"
list.fit[[4]]$VantHoff$Helix = "I"
list.fit[[4]]$VantHoff$Condition = "16mMMg-NTPs"

#NTPCM Helix J

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[5]] = read.csv("Figures/Figure_3/Fluorescence_data/js5057_Helix_J_formatted.csv")

list.fit[[5]] = meltR.F(list.df[[5]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/25mMNTPCM_Helix_J",
                        Save_results = "all",
                        K_error = c(0.62, 0.62),
                        B.conc.Tm = 200)

list.fit[[5]]$Tms$Helix = "J"
list.fit[[5]]$Tms$Condition = "25mMMg-NTPs"
list.fit[[5]]$VantHoff$Helix = "J"
list.fit[[5]]$VantHoff$Condition = "25mMMg-NTPs"

#2mMFree Helix F

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[6]] = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_F_formatted.csv")

list.fit[[6]] = meltR.F(list.df[[6]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/2mMFreeMg_Helix_F",
                        Save_results = "all",
                        B.conc.Tm = 200)

list.fit[[6]]$Tms$Helix = "F"
list.fit[[6]]$Tms$Condition = "2mMMg-free"
list.fit[[6]]$VantHoff$Helix = "F"
list.fit[[6]]$VantHoff$Condition = "2mMMg-free"

#2mMFree Helix G

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[7]] = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_G_formatted.csv")

list.fit[[7]] = meltR.F(list.df[[7]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/2mMFreeMg_Helix_G",
                        Save_results = "all",
                        B.conc.Tm = 300)

list.fit[[7]]$Tms$Helix = "G"
list.fit[[7]]$Tms$Condition = "2mMMg-free"
list.fit[[7]]$VantHoff$Helix = "G"
list.fit[[7]]$VantHoff$Condition = "2mMMg-free"

#2mMFree Helix H

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[8]] = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_H_formatted.csv")

list.fit[[8]] = meltR.F(list.df[[8]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/2mMFreeMg_Helix_H",
                        Save_results = "all",
                        B.conc.Tm = 300)

list.fit[[8]]$Tms$Helix = "H"
list.fit[[8]]$Tms$Condition = "2mMMg-free"
list.fit[[8]]$VantHoff$Helix = "H"
list.fit[[8]]$VantHoff$Condition = "2mMMg-free"


#2mMFree Helix I

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[9]] = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_I_formatted.csv")

list.fit[[9]] = meltR.F(list.df[[9]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/2mMFreeMg_Helix_I",
                        Save_results = "all",
                        B.conc.Tm = 300)

list.fit[[9]]$Tms$Helix = "I"
list.fit[[9]]$Tms$Condition = "2mMMg-free"
list.fit[[9]]$VantHoff$Helix = "I"
list.fit[[9]]$VantHoff$Condition = "2mMMg-free"

#2mMFree Helix J

paste("Figures/Figure_3/Fluorescence_data", list.files("Figures/Figure_3/Fluorescence_data"), sep = "/")

list.df[[10]] = read.csv("Figures/Figure_3/Fluorescence_data/js5060_Helix_J_formatted.csv")

list.fit[[10]] = meltR.F(list.df[[10]],
                        file_prefix = "Figures/Figure_3/Fluorescence_data/2mMFreeMg_Helix_J",
                        Save_results = "all",
                        B.conc.Tm = 300)

list.fit[[10]]$Tms$Helix = "J"
list.fit[[10]]$Tms$Condition = "2mMMg-free"
list.fit[[10]]$VantHoff$Helix = "J"
list.fit[[10]]$VantHoff$Condition = "2mMMg-free"

####Look at Tms####

df.Tms = bind_rows(list.fit[[1]]$Tms,
                   list.fit[[2]]$Tms,
                   list.fit[[3]]$Tms,
                   list.fit[[4]]$Tms,
                   list.fit[[5]]$Tms,
                   list.fit[[6]]$Tms,
                   list.fit[[7]]$Tms,
                   list.fit[[8]]$Tms,
                   list.fit[[9]]$Tms,
                   list.fit[[10]]$Tms)

head(df.Tms)

ggplot(df.Tms, aes(x = Condition, y = Tm, color = Helix)) +
  geom_beeswarm() +
  theme_classic()

####Look at free energies####


df.Energy = bind_rows(list.fit[[1]]$VantHoff,
                      list.fit[[2]]$VantHoff,
                      list.fit[[3]]$VantHoff,
                      list.fit[[4]]$VantHoff,
                      list.fit[[5]]$VantHoff,
                      list.fit[[6]]$VantHoff,
                      list.fit[[7]]$VantHoff,
                      list.fit[[8]]$VantHoff,
                      list.fit[[9]]$VantHoff,
                      list.fit[[10]]$VantHoff)

df.Energy

df.Energy = df.Energy %>%
  filter(Method == "2 Global fit") %>%
  arrange(Helix)

ggplot(df.Energy, aes(x = Helix, y = G, ymin = G - SE.G, ymax = G + SE.G, color = Condition))+
  geom_pointrange()

ggplot(df.Energy, aes(x = Helix, y = H, ymin = H - SE.H, ymax = H + SE.H, color = Condition))+
  geom_pointrange()

ggplot(df.Energy, aes(x = Helix, y = S, ymin = S - SE.S, ymax = S + SE.S, color = Condition))+
  geom_pointrange()
