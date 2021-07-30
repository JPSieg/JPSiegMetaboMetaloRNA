####Load the package####

devtools::load_all()

####Make a vector of data frames to read####

file.path.vector = paste("Data", list.files("Data"), sep = "/")

data.list = lapply(file.path.vector, read.itc)

head(data.list[[1]]$inj)
head(data.list[[5]]$inj)

?read.itc


####ATP 25 C####

df.M = read.itc("Data/M15mMMginto100uMATP25CpH7js2021.itc")
df.B = read.itc("Data/B15mMMgintoB25CpH7js2021.itc")


?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
            Save.path.and.prefix = "Results/M15mMMginto100uMATP25CpH7")

fit

####5 mM Mg into 0.5 mM EDTA####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/M5mMMginto500uMEDTA25CHEPESpH7js2024.itc", 5, 0, 0, 0, 0.5, 0))
df.B = read.itc("Data/B5mMMgintoB25CHEPESpH7js2024.itc")

?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/5_mM_Mg_into_0.5_mM_Mg",
                  Remove.injection = 1,
                  Fit.start = list(H = 306, K = 28800))

fit

####6 mM EDTA into 0.5 mM Mg ####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/M6mMEDTAinto500uMMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 0, 0.5, 0))
df.B = read.itc("Data/B6mMEDTAinto240mMNaCl140mMKCl10mMHEPES25CpH7.itc")

?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/6_mM_EDTA_into_0.5_mM_Mg",
                  Remove.injection = 1:5)

fit

####6 mM EDTA into Buffer####

?read.itc

read.itc("Data/B6mMEDTAinto240mMNaCl140mMKCl10mMHEPES25CpH7.itc", print.peak.integration.graph = TRUE, integration.folder = "Results")

####6 mM EDTA into 1.5 mM Mg 1 mM EDTA####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/M6mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 0))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")


?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/6_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA")

fit

####6 mM EDTA into 1.5 mM Mg 1 mM EDTA####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/M6mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 0))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")


?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/6_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA")

fit

####6 mM EDTA into 1.5 mM Mg 1 mM EDTA 150 mM Glutamate####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/C15mMEDTAinto500uMfreeMg150mMGlu1mMEDTAMg240mMNaCl140mMKCl25CpH7.itc", 15, 0, 0, 1, 1.5, 0))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")


?MetaboMgITC

fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/6_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA")

fit


####15 mM EDTA into 6 mM Mg 5 mM EDTA####

list.files("Data")

?read.itc

df.M = read.itc(file.concentrations = c("Data/M15mMEDTAinto6mMMg5mMEDTA240mNaCl140mMKCl25CpH7.itc", 14, 0, 0, 5, 6, 0))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")


fit = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/15_mM_EDTA_into_6_mM_Mg_5_mM_EDTA")

fit


####6 mM EDTA into 1.5 mM Mg 1 mM EDTA 150 mM Glutamate####


list.files("Data")

df.M = read.itc(file.concentrations = c("Data/M6mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 0))
df.C = read.itc(file.concentrations = c("Data/C6mMEDTAinto1.5mMMg1mMEDTA150mMGlutamate240mMNa140mMKCl10mMHEPESpH725C.itc", 6, 0, 150, 1, 1.5, 150))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")
df.C.B = read.itc(file.concentrations = c("B6mMEDTAinto25mMGlutathione1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 150, 1, 1.5, 150))

?MetaboMgITC.compete

fit = MetaboMgITC.compete(df.M, df.C, df.B, Save.path.and.prefix = "Results/Test")

fit



fit = MetaboMgITC(df.M, df.B)

fit

fit.C = MetaboMgITC(df.C, df.B,
                  Save.path.and.prefix = "Results/6_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA_150_mM_Glutamate")

fit.C

####6 mM EDTA into 1.5 mM Mg 1 mM EDTA 150 mM Glutamthione####


list.files("Data")

df.M = read.itc(file.concentrations = c("Data/M6mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 0))
df.C = read.itc(file.concentrations = c("Data/C6mMEDTAinto500uMMgFree25mMGlutathione1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 25))
df.B = read.itc("Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")
df.C.B = read.itc(file.concentrations = c("Data/B6mMEDTAinto25mMGlutathione1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 0, 1.5, 25))

?MetaboMgITC.compete

cell = df.M
competator = df.C
blank = df.B
blank.competator = df.C.B
Thermodynamic.equation = "Wiseman.isotherm"
Fit.start = list(V0 = 1.4247, H = 2534, K = 2880, Hc = 2534, Kc = 2880, c = 0, h = 0, Mg.contaminate = 1)
Remove.injection = 1
Saturation.threshold = FALSE
Save.path.and.prefix = FALSE

fit = MetaboMgITC.compete(df.M, df.C, df.B, Save.path.and.prefix = "Results/Test")

fit



fit = MetaboMgITC(df.M, df.B)

fit

fit.C = MetaboMgITC(df.C, df.B,
                    Save.path.and.prefix = "Results/6_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA_150_mM_Glutamate")

fit.C


####15 mM EDTA into 1.5 mM Mg 1 mM EDTA####

list.files("Data")



df.M = read.itc(file.concentrations = c("Data/M15mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl25CpH7.itc", 14, 0, 0, 1, 1.5, 0))
df.B = read.itc("Data/B15mMEDTAinto1mMEDTA240mMNaCl140mMKCl25CpH7.itc")

?MetaboMgITC

fit.M = MetaboMgITC(df.M, df.B,
                  Save.path.and.prefix = "Results/15_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA")

fit.M

####15 mM EDTA into 1.5 mM Mg 1 mM EDTA 150 mM Glutamate####

list.files("Data")

df.M = read.itc(file.concentrations = c("Data/C15mMEDTAinto500uMfreeMg150mMGlu1mMEDTAMg240mMNaCl140mMKCl25CpH7.itc", 15, 0, 0, 1, 1.5, 0))
df.B = read.itc("Data/B15mMEDTAinto150mMGlu1mMEDTA240mMNaCl140mMKCl25CpH7.itc")

?MetaboMgITC

fit.C = MetaboMgITC(df.M, df.B,
                    Save.path.and.prefix = "Results/15_mM_EDTA_into_1.5_mM_Mg_1_mM_EDTA")

fit.C


K = fit.M$Table[3,2]
K.app = fit.C$Table[3,2]

K.c = (K - K.app)/(K.app*0.15)

K.c

1000/0.81
