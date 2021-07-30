####Set up computer####

devtools::load_all()

library(tidyverse)

####Load in metabolite concentration data for E. coli####

list.files("Binding_constant_concentration_data")

df = read.csv("Binding_constant_concentration_data/210526_Metabolites.csv")

colnames(df)

unique(df$Species)

df.E.coli = df %>% filter(Species == "E.coli")

####Calculate Kd with high throughput####

?Kd.app.calc

View(df.E.coli)

IS = 0.390

Kd = c()

for (i in 1:length(df.E.coli$Metabolites)){
  Kd[i] = Kd.app.calc(df.E.coli$Metabolites[i], pH = 7.0, I = IS)
}

Kd

df.E.coli$Kd = Kd

colnames(df.E.coli)

####Write a csv to edit by hand####

?select

df.write = df.E.coli %>% select(Metabolites, Compartment, Concentration, Kd)

View(df.write)

write.csv(df.write, "E.coli_metabolites.csv", na = "", row.names = FALSE)
