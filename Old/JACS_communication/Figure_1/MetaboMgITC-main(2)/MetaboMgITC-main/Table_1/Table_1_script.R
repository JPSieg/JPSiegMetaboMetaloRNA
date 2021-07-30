devtools::document()
devtools::load_all()

?Kd.app.calc

####Load in Metabolite concentration data####

list.files("Binding_constant_concentration_data")

df = read.csv("Binding_constant_concentration_data/210526_Metabolites.csv")

df$Metabolites

####Calculate apparant Kds for IS = 0.15 and pH = 7.5####

IS = 0.15
pH = 7.5

#Test

Kd.app = Kd.app.calc("Lactic acid", pH, IS)
Kd.app

#Run in a for loop

Kd.app <- c()

for (i in 1:length(df$Metabolites)){
  print(i)
  print(df$Metabolites[i])
  Kd.app[i] = Kd.app.calc(df$Metabolites[i], pH, IS)
}

df$Kd.app = Kd.app
Source = rep(NA, length(Kd.app))
Source[which(is.na(Kd.app) == FALSE)] = "Calculated in this work"

df$Source = Source

####Print the table and edit by hand####

list.files()

write.csv(df, "Table 2/Table_2.csv", row.names = FALSE)

####Pull back in the hand curated csv and format it for Ryota####

list.files("Table 2")

df = read.csv("Table 2/Table_2_hand_edited.csv")

head(df)

df$MCM = df$Concentration*df$M/(df$M + df$Kd.app)

View(df)

df$Exp.biol.effect = df$MCM*df$Kd.app

colnames(df)

unique(df$Compartment)

df = dplyr::filter(df, Compartment != "Mitochondria")

unique(df$Compartment)

df = dplyr::select(df, Species, Metabolites, Concentration, Kd.app, MCM, Exp.biol.effect, Source, Notes)

write.csv(df, "Table_2_final.csv", row.names = FALSE)

head(df)
