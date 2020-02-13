# description: "Yeast specific transformation file (part of src/transform.R)"
# author: "Manuel Belmadani"
# date: "29/04/2018"
# updated: "29/04/2018"

yeast.tf.normalize_log2ratio <- 
  # Normalizes all variants based on empty vectors and wildtypes averages
  # For a spot Xi of a given variant:
  # where _Wt is the mean of WT and _Ev is the mean of Empty Vectors for a given sentinel
  # Compute the normalized score Xi``
  # Xi' = Xi / ( _Wt - _Ev )
  # Xi'' = Xi' - (_Wt - 1 )
  function(X, variable="log2ratio", STR_EMPTYVECTOR = "Empty vector") {
    
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    boundFunc <- function(Xi) {return (Xi - (right - 1)) }
    
    for(sent in unique(X$Gene.Name)){
      filt <- X$Gene.Name == sent
      left <- mean(na.rm=T, X$log2ratio[ (X$Class == STR_EMPTYVECTOR) & filt] )
      right <- mean( X$log2ratio[ (X$Class == STR_WT) & filt] )
      X$normalized[filt] <- sapply(X$log2ratio[filt], FUN = scaleFunc)
      X$normalized[filt] <- sapply(X$normalized[filt], FUN = boundFunc)
    }
    
    return(X)
  }

yeast.tf.normalizeByWT <- 
  # Scale WT to 1.0, leave control at 0.0
  function(X, variable, WT_CLASS = "WT", EV_CLASS = "EmptyVector") {
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    rownames(X) <- X$PTENVariant
    
    left <- 0
    right <- mean( X[ grepl(pattern=WT_CLASS, x=X$PTENVariant), variable] )
    
    X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
    return( X[, variable] )
}


yeast.tf.dropVariants <- 
  # Exclude variant based on known failed batches
  #   df - The entire variant dataset
  #   batch - The batch of the variants.
  #       1 - Original miniarray v2 (2015-2017)
  #       2 - Summer 2018 batch
  function(df, batch=1) {
    if (batch == 1){
      # Y180H, E285X and Q396R are excluded because the wrong construct was made.
      EXCLUSION_LIST = c("Y180H", "E285X", "Q396R")
      return(
        df[df$Variant %!in%  EXCLUSION_LIST,]
      )
      
    }
    if (batch == 2){
      # E307Q should be relabled as A309S 
      ## TODO: Make this part of fixTypos
      df[df$Variant == "E307Q","Variant"] <- c("A309S")
      # M35V and E157G did not mate and should be dropped
      EXCLUSION_LIST = c("M35V", "E157G")
      return(
        df[df$Variant %!in%  EXCLUSION_LIST,]
      )
      
    }
    
  }
