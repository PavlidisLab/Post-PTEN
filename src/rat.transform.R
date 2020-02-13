
rat.tf.normalizeByWT <- function(X, variable, WT_CLASS = "Wildtype overexpression", EV_CLASS = "EmptyVector") {
  scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
  boundFunc <- function(Xi) {return (Xi) - (right -1) }
  
  
  WT_mean <- mean( X[ (X$Class == WT_CLASS) , variable] )
  #right <- mean( X[ (X$Class == EV_CLASS), variable] )
  left = 0 
  right = WT_mean
  
  X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
  
  # # Sanity check:
  # stopifnot((((TRUE == cor(X$scaled, X$log2ratio)) == scaleFunc(right)) == (1 - scaleFunc(left))))
  
  return( X[, variable] )
  
}


fixRatDF <- function(R, normalize=T){
  
  
  row.names(R) <- gsub(row.names(R), pattern = "_OE", fixed = T, replacement = "")
  R$PTENVariant <- row.names(R)
  R <- merge(R, ClassInformation, all.x =T, all.y=F)
  R[R$PTENVariant == "Hum_WT", ]$Class <- "Wildtype overexpression"
  
  rownames(R) <- R$PTENVariant
  if (normalize) {
    R$PSD.95.Density = rat.tf.normalizeByWT(R, variable = "PSD.95.Density")
    R$Gephyrin.Density = rat.tf.normalizeByWT(R, variable = "Gephyrin.Density")
    R$Total.Dendrite.Length = rat.tf.normalizeByWT(R, variable = "Total.Dendrite.Length")
    R$Soma.Size = rat.tf.normalizeByWT(R, variable = "Soma.Size")
  }
  return(R)
}
