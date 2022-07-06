Radius = 50

N.mol = 0.195*(4/3)*(Radius^3)*(10^-27)*(6.022*10^23)  

N.Mg = 0.031*(4/3)*(Radius^3)*(10^-27)*(6.022*10^23)

N.Glu = ceiling(96/195*N.mol)  

N.Fru = floor(15/195*N.mol)

N.NTP = floor((9.63+8.29+4.87+4.62)/195*N.mol)

N.GluT = floor(16/195*N.mol)

sum(N.Glu, N.GluT, N.Fru, N.NTP) - ceiling(N.mol)
