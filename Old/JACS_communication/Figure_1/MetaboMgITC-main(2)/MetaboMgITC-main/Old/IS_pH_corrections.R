devtools::document()
devtools::load_all()

library(tidyverse)

####EDTA####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01)

?log10K.IS.calculator.2

#Debye-Huckel treatment from https://doi.org/10.1016/j.talanta.2007.06.035

log10.K.IS =  log10K.IS.calculator.2(8.79, 0.1, IS, 2, -4, -2)

log10.K.IS

pKa1 = 10.17 + 0.15

pKa2 = 6.11 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

a2 = 1/(1 + 10^(7-pKa1) + 10^(pKa2 - 7))

1000/(a1*10^log10.K.IS)

x = c(1:1500)
y = c()

for (i in x){
  y[i] = log10K.IS.calculator.2(8.79, 0.1, x[i]/1000, 2, -4, -2)
}

df1 = data.frame("Ionic.Strength" = x/1000,
                 "log10.K" = y,
                 "Method" = "Method 2")

#Debye-Huckel treatment from Buffers for pH and metal ion control

log10K.0 = log10K.0.calculator(8.79, 0.1, 2, -4, -2)

log10K.0

?log10K.IS.calculator

log10.K.IS =  log10K.IS.calculator(8.79, 0.1, 0.1, 2, -4, -2)

log10.K.IS

log10.K.IS =  log10K.IS.calculator(8.79, 0.1, IS, 2, -4, -2)

log10.K.IS

x = c(1:1500)
y = c()

for (i in x){
  y[i] = log10K.IS.calculator(8.79, 0.1, x[i]/1000, 2, -4, -2)
}

df2 = data.frame("Ionic.Strength" = x/1000,
                 "log10.K" = y,
                 "Method" = "Method 1")

df = bind_rows(df1, df2)

head(df2)

ggplot(df, aes(x = Ionic.Strength, y = log10.K, color = Method)) +
  geom_line() +
  xlab("Ionic strength (M)") +
  ylab("Log10[K (1/M)]") +
  geom_vline(xintercept = c(0.1, IS)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic()

ggsave("Analysis_of_Debye_Huckel_methods.png")

####ATP####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01)

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(4.3, 0.1, 0.1, 2, -4, -2)

log10.K.IS

log10.K.IS =  log10K.IS.calculator.2(4.3, 0.1, IS, 2, -4, -2)

log10.K.IS

pKa1 = 6.53 + 0.15

pKa2 = 4 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

K.app = a1*10^log10.K.IS

####Glu 6P####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01)

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(2.47, 0.1, IS, 2, -2, 0)

log10.K.IS

pKa1 = 6.5 + 0.15

pKa2 = 1.46 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

a2 = 1/(1 + 10^(7-pKa1) + 10^(pKa2 - 7))

a2

K.app = a1*10^log10.K.IS

1000/K.app

####Glutamate####

#Kd1

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01+0.1+0.1)

log10.K.IS =  log10K.IS.calculator.2(1.9, 0.1, IS, 2, -2, 0)

log10.K.IS

pKa1 = 9.59 + 0.15

pKa2 = 4.2 + 0.15

#pH 7

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

K.app = a1*10^log10.K.IS

1000/K.app

#pH 8

a1 = 1/(1 + 10^(pKa1 - 8) + 10^(pKa2 + pKa1 - 16))

a1

K.app = a1*10^log10.K.IS

1000/K.app

log10.K.IS =  log10K.IS.calculator.2(1.17, 0, IS, 2, -1, 1)

log10.K.IS

#K2

log10.K.IS =  log10K.IS.calculator.2(1.17, 0, IS, 2, -1, 1)

log10.K.IS

#pH 7

a2 = 1/(1 + 10^(7-pKa1) + 10^(pKa2 - 7))

a2

K.app = a2*10^log10.K.IS

1000/K.app

#pH 8

a2 = 1/(1 + 10^(8-pKa1) + 10^(pKa2 - 8))

a2

K.app = a2*10^log10.K.IS

1000/K.app

