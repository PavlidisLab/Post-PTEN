
FIELDS = 
  c(
    "PTENVariant",
    "Class",
    "Yeast.activity",
    "Fly.activity",
    "Worm.activity",
    "Axon.activity",
    "PSD.95.Density",
    "Total.Dendrite.Length",
    "Soma.Size",
    "Gephyrin.Density"
  )

reduce.data.frames <- function(df.list){
  Reduce(function(x, y) merge(x, y, all=TRUE), df.list)
}


fly.normalized.df$Class[is.na(fly.normalized.df$Class)] <- "Unclassified"
activity.list = list(
  yeast.normalized.df,
  fly.normalized.df,
  worm.normalized.df,
  rat.normalized.df,
  axon.normalized.df
)

activity.ci.list = list(
  yeast.ci.normalized.df,
  fly.ci.normalized.df,
  worm.ci.normalized.df,
  rat.ci.normalized.df,
  axon.ci.normalized.df
)

activity.err.list = list(
  yeast.err.normalized.df,
  fly.err.normalized.df,
  worm.err.normalized.df,
  rat.err.normalized.df,
  axon.err.normalized.df
)

activity.df = reduce.data.frames(activity.list)[,FIELDS]
activity.ci.df = reduce.data.frames(activity.ci.list)[,FIELDS]
activity.err.df = reduce.data.frames(activity.err.list)[,FIELDS]

dim(activity.df)
dim(activity.ci.df)
dim(activity.err.df)
dim(protein.df)
dim(protein.err)

write.csv(x = activity.df, file = "Distributable/activity.csv", row.names = F)
write.csv(x = activity.ci.df, file = "Distributable/activity.ci.csv", row.names = F)
write.csv(x = activity.err.df, file = "Distributable/activity.err.csv", row.names = F)

write.csv(x = protein.df, file = "Distributable/protein.csv", row.names = F)
write.csv(x = protein.err, file = "Distributable/protein.err.csv", row.names = F)
