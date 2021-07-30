devtools::document()
devtools::load_all()

library(tidyverse)
library(cowplot)

####ITC files####

list.files("ITC_Data")

####EDTA####

#140 K pH 7.5

df.M.EDTA.140K.pH.7.5 = read.itc(c("ITC_Data/M5mMEDTAinto500uMfreeMg1mMEDTAMG140mMKCl10pMHEPESpH7p5.itc", 5, 0, 0, 1, 1.5, 0))
df.B.EDTA.140K.pH.7.5 = read.itc("ITC_Data/B5mMEDTAinto1mMEDTA140mMKCl10pMHEPESpH7p5.itc")

fit.EDTA.140K.pH.7.5 = MetaboMgITC(df.M.EDTA.140K.pH.7.5,
                                   df.B.EDTA.140K.pH.7.5,
                                   Fit.start = list(H = 3.6, K = 28800))
#240 Na 140 K pH 7.0

list.files("ITC_Data")

df.M.EDTA.240Na.140K.pH.7 = read.itc(c("ITC_Data/M6mMEDTAinto500uMfreeMg1mMEDTAMg240mMNaCl140mMKCl10mMHEPES25CpH7.itc", 6, 0, 0, 1, 1.5, 0))
df.B.EDTA.240Na.140K.pH.7 = read.itc("ITC_Data/B6mMEDTAinto1mMEDTA240mMNaCl140mMKCl10mMHEPES25CpH7.itc")

fit.EDTA.240Na.140K.pH.7 = MetaboMgITC(df.M.EDTA.240Na.140K.pH.7,
                                   df.B.EDTA.240Na.140K.pH.7,
                                   Fit.start = list(H = 3.6, K = 194))

####ATP####

#140 K pH 7.5

list.files("ITC_Data")

df.M.ATP.140K.pH.7.5 = read.itc(c("ITC_Data/M15mMMginto200uMATP140mMKCl10mMHEPESpH7p5.itc", 15, 0, 0, 0, 0.2, 0))
df.B.ATP.140K.pH.7.5 = read.itc("ITC_Data/B15mMMginto140mMKCl10mMHEPESpH7p5.itc")

fit.ATP.140K.pH.7.5 = MetaboMgITC(df.M.ATP.140K.pH.7.5,
                                   df.B.ATP.140K.pH.7.5,
                                   Fit.start = list(H = 3.6, K = 28800))
#240 Na 140 K pH 7.0

list.files("ITC_Data")

df.M.ATP.240Na.140K.pH.7 = read.itc(c("ITC_Data/M15mMMginto100uMATP25CpH7js2021.itc", 15, 0, 0, 0, 0.1, 0))
df.B.ATP.240Na.140K.pH.7 = read.itc("ITC_Data/B15mMMgintoB25CpH7js2021.itc")

fit.ATP.240Na.140K.pH.7 = MetaboMgITC(df.M.ATP.240Na.140K.pH.7,
                                   df.B.ATP.240Na.140K.pH.7,
                                   Fit.start = list(H = 3.6, K = 2880))

####Glucose 6-P####

#140 K pH 7.5

list.files("ITC_Data")

df.M.Glu6P.140K.pH.7.5 = read.itc(c("ITC_Data/M100mMMginto5mMGlu6P140mMKCl10mMHEPESpH7p5.itc", 100, 0, 0, 0, 5, 0))
df.B.Glu6P.140K.pH.7.5 = read.itc("ITC_Data/B100mMMgintoB25CHEPESpH7js2024.itc")

fit.Glu6P.140K.pH.7.5 = MetaboMgITC(df.M.Glu6P.140K.pH.7.5,
                                  df.B.Glu6P.140K.pH.7.5,
                                  Fit.start = list(H = 1, K = 5))

#240 Na 140 K pH 7.0

list.files("ITC_Data")

df.M.Glu6P.240Na.140K.pH.7 = read.itc(c("ITC_Data/M100mMMginto5mMGlu6P25CHEPESpH7js2024.itc", 100, 0, 0, 0, 5, 0))
df.B.Glu6P.240Na.140K.pH.7 = read.itc("ITC_Data/B100mMMgintoB25CHEPESpH7js2024.itc")

fit.Glu6P.240Na.140K.pH.7 = MetaboMgITC(df.M.Glu6P.240Na.140K.pH.7,
                                df.B.Glu6P.240Na.140K.pH.7,
                                Fit.start = list(H = 0.1, K = 5))

####Consolidate results####

#SI figure X

SI_figure_X = plot_grid(fit.EDTA.140K.pH.7.5$Plot, fit.EDTA.240Na.140K.pH.7$Plot)

list.files("SI_ITC_figures")

ggsave("SI_ITC_figures/SI_figure_1_EDTA.svg", SI_figure_X)

#SI figure X plus 1

SI_figure_X_plus_1 = plot_grid(fit.ATP.140K.pH.7.5$Plot, fit.ATP.240Na.140K.pH.7$Plot)

list.files("SI_ITC_figures")

ggsave("SI_ITC_figures/SI_figure_2_ATP.svg", SI_figure_X_plus_1)

#SI figure Y

SI_figure_Y = plot_grid(fit.Glu6P.140K.pH.7.5$Plot, fit.Glu6P.240Na.140K.pH.7$Plot)

list.files("SI_ITC_figures")

ggsave("SI_ITC_figures/SI_figure_3_Glucose_6P.svg", SI_figure_Y)
