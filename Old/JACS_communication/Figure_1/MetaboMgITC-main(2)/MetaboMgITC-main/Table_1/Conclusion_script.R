list.files("Table 1")

library(tidyverse)

df = read.csv("Table 1/Table_1_final.csv")

head(df)

unique(df$Species)

df = subset(df, df$Species == "E.coli")

sum(df$Concentration)

sum(df$MCM, na.rm = TRUE)

colnames(df)

df = df %>% filter(Kd.app >= 2)

sum(df$MCM, na.rm = TRUE)

df = read.csv("Table 1/Table_1_final.csv")

unique(df$Species)

df = subset(df, df$Species == "Mammal iBMK" )

sum(df$Concentration)

sum(df$MCM, na.rm = TRUE)

colnames(df)

df = df %>% filter(Kd.app >= 2)

sum(df$MCM, na.rm = TRUE)

df = read.csv("Table 1/Table_1_final.csv")

unique(df$Species)

df = subset(df, df$Species == "Yeast")

sum(df$Concentration)

sum(df$MCM, na.rm = TRUE)

colnames(df)

df = df %>% filter(Kd.app >= 2)

sum(df$MCM, na.rm = TRUE)
