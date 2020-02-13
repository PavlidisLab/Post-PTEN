# description: "Worm specific transformation file (part of src/transform.R)"
# author: "Manuel Belmadani"
# date: "29/04/2018"
# updated: "29/04/2018"

fly.tf.normalizeByWT  <- 
  # Descriptions
  function(X, variable, WT_VAR = "WT", EV_VAR = "attp2") {
    ## X=worm.df; variable="Worm.activity"
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    boundFunc <- function(Xi) {return ( (right - left) - Xi ) }
    
    rownames(X) <- X$PTENVariant
    
    left <- mean( X[ EV_VAR, variable] )
    right <- mean( X[ WT_VAR, variable] )
    
    X[, variable] = X[, variable] - left
    X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
    
    return( X[, variable] )
  }
