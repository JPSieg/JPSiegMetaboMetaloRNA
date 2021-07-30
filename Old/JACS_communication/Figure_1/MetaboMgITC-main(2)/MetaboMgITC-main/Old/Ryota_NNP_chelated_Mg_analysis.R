devtools::document()
devtools::load_all()

library(tidyverse)

####ATP####

IS = 0.165

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(4.3, 0.1, 0.1, 2, -4, -2)

log10.K.IS

log10.K.IS =  log10K.IS.calculator.2(4.3, 0.1, IS, 2, -4, -2)

log10.K.IS

pKa1 = 6.53 + 0.15

pKa2 = 4 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7.5) + 10^(pKa2 + pKa1 - 15))

a1

K.app = a1*10^log10.K.IS

1000/K.app


X.0 = 9.6
M.0 = 8.7
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM


X.0 = 9.6
M.0 = 11.2
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

####ADP####

IS = 0.165

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(3.17, 0.1, IS, 2, -3, -1)

log10.K.IS

pKa1 = 6.4 + 0.15

pKa2 = 3.96 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7.5) + 10^(pKa2 + pKa1 - 15))

a1

K.app = a1*10^log10.K.IS

1000/K.app

X.0 = 0.56
M.0 = 0.74
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM


X.0 = 0.56
M.0 = 2.42
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

####AMP####

IS = 0.165

?log10K.0.calculator

log10.K.IS =  log10K.IS.calculator.2(1.93, 0.1, IS, 2, -2, 0)

log10.K.IS

pKa1 = 6.02 + 0.15

pKa2 = 3.72 + 0.15

a1 = 1/(1 + 10^(pKa1 - 7.5) + 10^(pKa2 + pKa1 - 15))

a1

K.app = a1*10^log10.K.IS

1000/K.app

X.0 = 0.28
M.0 = 0.501
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM


X.0 = 0.28
M.0 = 2.04
Kd = 1000/K.app

b = -(X.0 + Kd + M.0)
c = X.0*M.0
XM = (-b - sqrt((b)^2 - 4*c))/2
XM

