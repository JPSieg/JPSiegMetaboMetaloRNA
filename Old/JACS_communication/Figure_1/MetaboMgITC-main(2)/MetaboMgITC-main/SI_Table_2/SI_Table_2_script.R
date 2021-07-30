devtools::document()
devtools::load_all()

library(tidyverse)

####Theoretical K app####

Kd.app.EDTA.0.15.pH.7.5 = Kd.app.calc("EDTA", I = 0.15)
Kd.app.EDTA.0.15.pH.7.5
Kd.app.EDTA.0.39.pH.7.0 = Kd.app.calc("EDTA", pH = 7.0, I = 0.39)
Kd.app.EDTA.0.39.pH.7.0
Kd.app.EDTA.0.15.pH.12 = Kd.app.calc("EDTA", pH = 12, I = 0.15)
Kd.app.EDTA.0.15.pH.12
Kd.app.EDTA.0.39.pH.12 = Kd.app.calc("EDTA", pH = 12, I = 0.39)
Kd.app.EDTA.0.39.pH.12


Kd.app.ATP.0.15.pH.7.5 = Kd.app.calc("ATP", pH = 7.5, I = 0.15)
Kd.app.ATP.0.15.pH.7.5
Kd.app.ATP.0.39.pH.7.0 = Kd.app.calc("ATP", pH = 7.0, I = 0.39)
Kd.app.ATP.0.39.pH.7.0
Kd.app.ATP.0.15.pH.12 = Kd.app.calc("ATP", pH = 12, I = 0.15)
Kd.app.ATP.0.15.pH.12
Kd.app.ATP.0.39.pH.12 = Kd.app.calc("ATP", pH = 12, I = 0.39)
Kd.app.ATP.0.39.pH.12

Kd.app.Glucose.6P.0.15.pH.7.5 = Kd.app.calc("Glucose 6-P", I = 0.15)
Kd.app.Glucose.6P.0.15.pH.7.5
Kd.app.Glucose.6P.0.39.pH.7.0 = Kd.app.calc("Glucose 6-P", pH = 7.0, I = 0.39)
Kd.app.Glucose.6P.0.39.pH.7.0
Kd.app.Glucose.6P.0.15.pH.12 = Kd.app.calc("Glucose 6-P", pH = 12, I = 0.15)
Kd.app.Glucose.6P.0.15.pH.12
Kd.app.Glucose.6P.0.39.pH.12 = Kd.app.calc("Glucose 6-P", pH = 12, I = 0.39)
Kd.app.Glucose.6P.0.39.pH.12



####ITC files####

list.files("ITC_Data")

?read.itc

####EDTA####

#140 K pH 7.5

df.M.EDTA.140K.pH.7.5 = read.itc(c("ITC_Data/M5mMEDTAinto500uMfreeMg1mMEDTAMG140mMKCl10pMHEPESpH7p5.itc", 5, 0, 0, 1, 1.5, 0))
df.B.EDTA.140K.pH.7.5 = read.itc("ITC_Data/B5mMEDTAinto1mMEDTA140mMKCl10pMHEPESpH7p5.itc")

?MetaboMgITC

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

?MetaboMgITC

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

?MetaboMgITC

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

(0.14 + 0.01)
(0.240 + 0.14 + 0.01)


Metabolite = c("EDTA", "EDTA", "ATP", "ATP", "Glucose 6-P", "Glucose 6-P")
IS = c(0.15, 0.39, 0.15, 0.39, 0.15, 0.39)
pH = c(7.5, 7.0, 7.5, 7.0, 7.5, 7.0)
ITC.Kd.app = c(1000/fit.EDTA.140K.pH.7.5$Table$Number[3],
               1000/fit.EDTA.240Na.140K.pH.7$Table$Number[3],
               1000/fit.ATP.140K.pH.7.5$Table$Number[3],
               1000/fit.ATP.240Na.140K.pH.7$Table$Number[3],
               1000/fit.Glu6P.140K.pH.7.5$Table$Number[3],
               1000/fit.Glu6P.240Na.140K.pH.7$Table$Number[3])
Theoretical.Kd.app = c(Kd.app.EDTA.0.15.pH.7.5,
                       Kd.app.EDTA.0.39.pH.7.0,
                       Kd.app.ATP.0.15.pH.7.5,
                       Kd.app.ATP.0.39.pH.7.0,
                       Kd.app.Glucose.6P.0.15.pH.7.5,
                       Kd.app.Glucose.6P.0.39.pH.7.0)

High.pH.Kapp = c(Kd.app.EDTA.0.15.pH.12,
                 Kd.app.EDTA.0.39.pH.12,
                 Kd.app.ATP.0.15.pH.12,
                 Kd.app.ATP.0.39.pH.12,
                 Kd.app.Glucose.6P.0.15.pH.12,
                 Kd.app.Glucose.6P.0.39.pH.12)

No.IS.or.pK.corr = c(1000/(10^8.79),
                     1000/(10^8.79),
                     1000/(10^4.3),
                     1000/(10^4.3),
                     1000/(10^2.47),
                     1000/(10^2.47))


df = data.frame(Metabolite,
                No.IS.or.pK.corr,
                High.pH.Kapp,
                IS,
                pH,
                Theoretical.Kd.app,
                ITC.Kd.app)

df


write.csv(df, "SI_Table_2/SI_Table_2.csv", row.names = FALSE)
