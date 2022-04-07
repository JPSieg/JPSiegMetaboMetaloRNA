library(R2easyR)

####Gua aptamer####

vector.files = paste("Figures/Figure_2/Gua_data_files", list.files("Figures/Figure_2/Gua_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)
df = bind_rows(list.df)
hist(df$Reactivity, breaks =20)
abline(v = c(100, 1000))

palettes = r2easyR.palettes()
Condition = c("25mMFree", "2mMFree", "Eco80", "NTPCM", "WMCM")

table(list.df[[1]]$Dotbracket)

for (i in 1:length(list.df)){
  list.df[[i]]$Reactivity[list.df[[i]]$Reactivity <= 0] = NA
  list.df[[i]] = r2easyR.color(list.df[[i]],
                               palettes$Reds.c,
                               manual.scale = c(100, 1000))

  prefix = gsub("Data_files", "R2R", strsplit(vector.files[i], split = "[.]")[[1]][1])

  r2easyR.write(prefix, list.df[[i]], RNA_name = Condition[i], colors = "circles")
}

r2easyR.write("Figures/Figure_2/R2R/No_reactivity", list.df[[1]], RNA_name = "No_reactivity", colors = "NA")

####CPEB3 aptamer####

vector.files = paste("Figures/Figure_2/CPEB3_data_files", list.files("Figures/Figure_2/CPEB3_data_files"), sep = "/")

list.df = lapply(vector.files, read.csv)
df = bind_rows(list.df)
hist(df$Reactivity, breaks =20)
abline(v = c(100, 1000))

palettes = r2easyR.palettes()
Condition = c("25mMFree", "2mMFree", "Eco80", "NTPCM", "WMCM")

table(list.df[[1]]$Dotbracket)

for (i in 1:length(list.df)){
  list.df[[i]]$Reactivity[list.df[[i]]$Reactivity <= 0] = NA
  list.df[[i]] = r2easyR.color(list.df[[i]],
                               palettes$Reds.c,
                               manual.scale = c(100, 1000))

  pk1 = which(list.df[[i]]$Pknot == 1)
  pk3 = which(list.df[[i]]$Pknot == 3)

  list.df[[i]]$Dotbracket[c(pk1, pk3)] = "."

  prefix = paste(gsub("CPEB3_data_files", "R2R", strsplit(vector.files[i], split = "[.]")[[1]][1]), "_CPEB3", sep = "")

  r2easyR.write(prefix, list.df[[i]], RNA_name = Condition[i], colors = "circles")
  r2easyR.stem_editor(paste(prefix, ".sto", sep = ""))
  r2easyR.pknot_drawer(paste(prefix, ".sto", sep = ""), list(pk1, pk2, pk3))
}

r2easyR.write("Figures/Figure_2/R2R/CPEB3_no_reactivity", list.df[[1]], RNA_name = "No_reactivity", colors = "NA")
r2easyR.stem_editor("Figures/Figure_2/R2R/CPEB3_no_reactivity.sto")
r2easyR.pknot_drawer("Figures/Figure_2/R2R/CPEB3_no_reactivity.sto", list(pk1, pk3))
