library(tidyverse)

df.met = read.csv("Figures/Figure_1/Top_15_E.coli_metabolites_edited.csv")

models = 1000

AC.met.model = function(data_frame = df.met){
  M.error = c()
  Kd.error = c()
  for ( i in 1:nrow(data_frame)){
    M.error[i] = rnorm(1, 0, data_frame$M.error[i])
    #M.error = M.error + rnorm(1, 0.015, 0.005)
    if (is.na(data_frame$Kd.error[i])){}else{
      Kd.error[i] = rnorm(1, 0, data_frame$Kd.error[i])
    }
  }
  data_frame$Model.conc = data_frame$Concentration + M.error
  data_frame$Model.Kd = data_frame$Kd + Kd.error
  output = data_frame  %>% select(Metabolites,  Model.conc, Model.Kd, Mg.binding.strength)
}

Numeric.free.Mg.solver = function(df, Mg.total = 4.042252){
  go = TRUE
  First.recursion = TRUE

  while(go){
    if (First.recursion){
      Mg.free = seq(0.000000001, Mg.total, length.out = 20)
      First.recursion = FALSE
    }

    error.Mg.T = c()

    for (i in 1:length(Mg.free)){
      MCM = df$Model.conc*Mg.free[i]/(df$Model.Kd + Mg.free[i])
      error.Mg.T[i] = abs(Mg.total - Mg.free[i] - sum(MCM, na.rm = TRUE))
    }

    center = Mg.free[which.min(error.Mg.T)]

    error = error.Mg.T[which.min(error.Mg.T)]

    if (center - error <= 0){
      Mg.free = seq(0, center + error, length.out = 100)
    }else{
      Mg.free = seq(center - error, center + error, length.out = 50)
    }

    if (error <= 0.01){go = FALSE}
  }
  #print(center)
  output = center
}

list.df.model = {}
Iteration = c()
MCM = c()
MCM.conc = c()

pb = txtProgressBar(min = 1, max = models, initial = 1)

print("Modeling artificial cytoplasms")

for (i in 1:models){
  setTxtProgressBar(pb, i)
  list.df.model[[i]] = AC.met.model()
  list.df.model[[i]]$Iteration = i
  df.NTPCM = list.df.model[[i]] %>% filter(Mg.binding.strength == "strong")
  Iteration = c(Iteration, i)
  MCM = c(MCM, "NTPCM")
  MCM.conc = c(MCM.conc, sum(df.NTPCM$Model.conc))
  df.WMCM = list.df.model[[i]] %>% filter(Mg.binding.strength == "weak")
  Iteration = c(Iteration, i)
  MCM = c(MCM, "WMCM")
  MCM.conc = c(MCM.conc, sum(df.WMCM$Model.conc))
  Iteration = c(Iteration, i)
  MCM = c(MCM, "Ecoli80")
  MCM.conc = c(MCM.conc, sum(list.df.model[[i]]$Model.conc))
  close(pb)
}


df = data.frame(Iteration,
                MCM,
                MCM.conc)

write.csv(df, "Figures/Figure_1/Modeled_AC_metabolite_concentrations.csv", row.names = FALSE)



Mg.T = 10^seq(-1, log10(200), length.out = 250)

pb = txtProgressBar(min = 1, max = length(list.df.model), initial = 1)

list.df.MCM.model = {}

print("Modeling free Mg concentrations")

for (i in 1:length(list.df.model)){
  setTxtProgressBar(pb, i)
  Mg.free =  c()
  Mg.free.NTP = c()
  Mg.free.WMCM = c()
  for (j in 1:length(Mg.T)){
    #print(j)
    #Total MCM
    df = list.df.model[[i]]
    Mg.free[j] = Numeric.free.Mg.solver(df, Mg.T[j])
    #Total MCM
    df = list.df.model[[i]] %>% filter(Mg.binding.strength == "strong")
    Mg.free.NTP[j] = Numeric.free.Mg.solver(df, Mg.T[j])
    #Total MCM
    df = list.df.model[[i]] %>% filter(Mg.binding.strength == "weak")
    Mg.free.WMCM[j] = Numeric.free.Mg.solver(df, Mg.T[j])
  }

  list.df.MCM.model[[i]] = data.frame(i,
                                      Mg.T,
                                      Mg.free,
                                      Mg.free.NTP,
                                      Mg.free.WMCM)

  close(pb)
}

df = bind_rows(list.df.MCM.model)
write.csv(df, "Figures/Figure_1/Modeled_AC_MCM_concentrations.csv", row.names = FALSE)

