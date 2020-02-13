rat.raw = read.csv(sep="\t",
              "RawData/Rat/2018-11-05_PTEN-Rat.tsv")[-1,1:7] # Subset for data rows only


rat.raw$Condition = gsub(pattern = " ", replacement = "_", rat.raw$Condition)
rat.raw$Condition = gsub(pattern = ".*\\._", replacement = "", rat.raw$Condition)
rat.raw$Condition = gsub(pattern = "_OE", replacement = "", rat.raw$Condition)

rat.raw$Cell.ID = gsub(pattern = " ", replacement = "_", rat.raw$Cell.ID)
rat.raw$Culture = gsub(pattern = " ", replacement = "_", rat.raw$Culture)

write.table(rat.raw, file = "Data/PTEN-Rat.tsv", sep = "\t")
