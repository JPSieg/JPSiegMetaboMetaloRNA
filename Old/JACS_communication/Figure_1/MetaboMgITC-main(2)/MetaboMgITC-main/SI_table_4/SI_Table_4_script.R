devtools::load_all()

Kd.app.calc("L-Glutamic acid",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("L-Aspartic acid",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("L-Glutamine",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("L-Alanine",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("ATP",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("ADP",
            pH = 7.5,
            I = 0.15)

Kd.app.calc("AMP",
            pH = 7.5,
            I = 0.15)

list.files("SI_table_4")

df = read.csv("SI_table_4/SI_Table_4_start.csv")

a = 1
b = df$L.total + df$M.total + df$Kd.IS.pH.corr
c = df$L.total*df$M.total

MCM.corr = (b - sqrt(b^2 - 4*a*c))/(2*a)

MCM.corr

M.free.corr = df$M.total - MCM.corr

M.free.corr

M.free.difference = M.free.corr - df$M.free.mM.original

M.free.difference

M.free.increase = M.free.corr/df$M.free.mM.original

M.free.increase

df$MCM.corr = MCM.corr
df$M.free.corr = M.free.corr
df$M.free.increase = M.free.increase

df

write.csv(df, "SI_table_4/SI_table_4.csv", row.names = FALSE)
