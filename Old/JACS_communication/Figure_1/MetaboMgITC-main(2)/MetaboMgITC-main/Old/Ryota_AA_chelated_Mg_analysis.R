devtools::document()
devtools::load_all()

library(tidyverse)

####EDTA####

IS = 0.165

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(8.79, 0.1, IS, 2, -4, -2)

log10.K.IS

pKa1 = 10.17 + 0.15

pKa2 = 6.11 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7.5) + 10^(pKa2 + pKa1 - 15))

a1

K.app = a1*10^log10.K.IS

1000/K.app

#0.1 mM Free Mg

X.0 = 0.66
M.0 = 0.76
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

#0.5 mM Free Mg

X.0 = 3.1
M.0 = 3.6
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

#2 mM Free Mg

X.0 = 11.3
M.0 = 13.3
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM


####Glu####

IS = 0.165

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(1.17, 0.1, IS, 2, -1, 1)

log10.K.IS

pKa1 = 9.59 + 0.15

pKa2 = 4.2 + 0.15

a2 = 1/(1 + 10^(7.5-pKa1) + 10^(pKa2 - 7.5))

a2

K.app = a2*10^log10.K.IS

1000/K.app

#0.1 mM Free Mg

X.0 = 96
M.0 = 0.76
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

#0.5 mM Free Mg

X.0 = 96
M.0 = 3.6
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

#2 mM Free Mg

X.0 = 96
M.0 = 13.3
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM
