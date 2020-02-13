# description: "Preprocessing the Raw PTEN yeast miniarray data; keep each data point separate but pair them with wildtypes."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
set.seed(42)
source("utils.R")

# SET CONFIGS
dataFilesPath.batch1 <- "RawData/YeastActivity/miniarray_v2"
dataFilesPath.batch2 <- "RawData/YeastActivity/180603_PTENMini_MB"

## Load data from per-plate csv
df.batch1 <- yeast.io.load_multiple_csv(dataFilesPath.batch1)
df.batch2 <- yeast.io.load_multiple_csv(dataFilesPath.batch2)

# Fixes
table(df.batch1$Variant)
table(df.batch2$Variant)

# Fix known typos
df.batch1$Variant <- yeast.io.fixTypos(x = df.batch1$Variant, batch=1)
df.batch2$Variant <- yeast.io.fixTypos(x = df.batch2$Variant, batch=2)

df.batch1 <- yeast.tf.dropVariants(df = df.batch1, batch=1)
df.batch2 <- yeast.tf.dropVariants(df = df.batch2, batch=2)

# Merge batches
df <- rbind(df.batch1, df.batch2)
df <- df[order(df$LogicalRow), ]

df <- yeast.io.match_pten_metadata(df, PTEN_METADATA) #
rownames(df) <- 1:nrow(df)

classesVariants<-unique(df[,c("Variant", "Class")])
classes<-classesVariants[,-1,drop=F]

classColorsOrder <- as.numeric(c())
for (name in names(classColors)) {
  # print(which( (name == classes$Class) == TRUE ))
  rows <- which( (name == classes$Class) == TRUE )
  classColorsOrder <- c( classColorsOrder, rows)
}

unique.spots = unique(df[,c( "Row", "Column", "block", "Screen", "Gene.Name", "Class", "Variant")])
df$paired <- NA
for ( i in 1:nrow(unique.spots)){
  indices = !is.na(prodlim::row.match(df[,names(unique.spots)], unique.spots[i,]))
  df[indices,]$paired <- i
}

head(df)

write.table(df,
            file = "Data/PTEN-MEM-Yeast.tsv", sep = '\t', row.names = F)
