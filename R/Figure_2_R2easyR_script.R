library(R2easyR)

vector.files = paste("Figures/Figure_2/Data_files", list.files("Figures/Figure_2/Data_files"), sep = "/")



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
