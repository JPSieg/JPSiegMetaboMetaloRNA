library(tidyverse)

list.files("Binding_constant_concentration_data")

df = read.csv("Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv")

head(df)

df = df[-which(is.na(df$Log_K)),]

head(df)

df = df %>% arrange(Metabolite)

View(df)

table(df$Citation)

pattern = unique(df$Citation)

replacement = c("[6]", "[9]", "[5]", "", "[4]", "[7]", "[10]", "[11]", "[12]", "[13]")

for (i in 1:length(df$Citation)){
  df$Citation[i] = replacement[which(pattern == df$Citation[i])]
}

View(df)

list.files("SI_Table_3")

write.csv(df, "SI_Table_3/SI_Table_3_final.csv", row.names = FALSE)
